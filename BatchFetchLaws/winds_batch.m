function [sigH,htgrid,E_each] = winds_batch(planet,liquid,model,wind,Etc,ann,bf1,bf2,c,ccw,Cg,cgmax,cm,cp,cth,cth2,cw,D,delx,dely,dr,dth,dthd,dwn,E,E_each,f,file,freqs,gust,htgrid,i,idx,kappa,kutoff,l2,lz,mindelx,mss_fac,nnn,nnninv,oa,ol,os,rhoa,rhorat,Sbf_fac,Sdt_fac,sigH,sm,Snl_fac,sp,sth,TitanResults,waveang,wfac,wn,xm,xp,ym,yp,z)

   UU = wind.speed;                                                        % UU = wind speed for loop
   file = file + 1;                                                        % for naming files
 
   U = UU*ones(model.m,model.n);                                           % set wind velocity everywhere in x-y plane
   windir = wind.dir*ones(model.m,model.n);                                % set wind direction everywhere in x-y plane
   windir = repmat(windir,[1 1 model.o model.p]);                          % reshape wind direction matrix
    
   U_z = U + 0.005;                                                        % wind at modelled height (plus small number to wind speed to avoid division by zero)
  
   % drag coefficient 
   Cd = 1.2*ones(size(U));                                                 % drag coefficient for weak/moderate winds
   Cd(U > 11) = 0.49 + 0.065*U(U > 11);                                    % drag coefficient for strong winds
   Cd = Cd/1000;
  
   modt = 0;                                                               % model time initializiation
  
   % Currents superimposed onto wave-driven flow
   Uer = 0*D + 0.0;                                                        % Eastward current [m/s]
   Uei = 0*D + 0.0;                                                        % Northward current [m/s]
  
%% -- loop through model time to grow waves --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   for t = 1:model.num_time_steps                                                                                                                                    % loop through time
      
       sumt = 0;                                                                                                                                                     % intitalize total time within timestep t
       tplot = - 1; 
       while ((model.time_step - sumt) > 0)                                                                                                                          % each time step is determined as min[max(0.1/(Sin-Sds) 0.0001) 2000 UserDefinedTime CourantGrid], iterate until the user-defined time is larger than the time passed within timestep t
          
           if gust > 0
                U = U.*(1+gust*randn(size(U)));                                                                                                                      % add random bursts of gusts to wind speed
           end

           tplot = tplot + 1;
           explim = 0.1;                                                                                                                                             % limit for making dynamic time step if the source = dissipation
           delt = 0.7 * mindelx / cgmax;                                                                                                                             % advection limited Courant condition.
          
           Ul = real(U_z.*sqrt(Cd));                                        
           Ustar = Ul;
           U_10 = 2.5*Ul.*log(10/z) + U_z;                                                                                                                           % law of wall
           ustarw = Ustar .* sqrt(rhorat);                                                                                                                           % shear stress
          
           Ul = repmat(Ul,[1 1 model.o model.p]);
           U = repmat(U_z,[1 1 model.o model.p]);
           Ul = 2.5*Ul.*log(l2/z)+U;

           ustw = repmat(ustarw,[1 1 model.o model.p]);
           Ud = - ustw./kappa.*log(lz/z);                                                                                                                            % depth averaged air velocity (law of the wall)
          
           % calculate wind input as a fraction of stress
           Ul(D>0) = abs(Ul(D>0).*cos(windir(D>0)-waveang(D>0))-c(D>0)-Ud(D>0).*cos(windir(D>0)-waveang(D>0))-Uer(D>0).*cth(D>0)-Uei(D>0).*sth(D>0))...              % Sin = A1*Ul*((k*wn)/g)*(rhoa/rhow)*E  [eqn. 4, Donelan 2012]
               .*(Ul(D>0).*cos(windir(D>0)-waveang(D>0))-c(D>0)-Ud(D>0).*cos(windir(D>0)-waveang(D>0))-Uer(D>0).*cth(D>0)-Uei(D>0).*sth(D>0));
          
          
           % Adjust input to lower value when waves overrun the wind because as wave speed approaches wind speed, energy extraction becomes less efficient
           Ula = Ul;
           Ul = 0.11*Ul;
           Ul(Ul<0) = Ula(Ul<0)*0.03;                                                                                                                                % Field scale (Donelan et al, 2012). 0.135 in Lab                        
           Ul(cos(windir-waveang)<0) = Ula(cos(windir-waveang)<0)*0.015;                                                                                             % Waves against wind

% -- Sin --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           % input to energy spectrum from wind
           Sin = zeros(size(Ul));
           Sin(D>0) = rhorat*(Ul(D>0).*2*pi.*f(D>0).*wn(D>0)./(planet.gravity+liquid.surface_tension.*wn(D>0).*wn(D>0)./liquid.rho_liquid));         % eqn. 4, Donelan 2012
           % limits energy going into spectrum once waves outrun wind
           Heavyside = ones(size(E));
           Heavyside(Sin < 0 & E < 1e-320) = 0;
           Sin = Heavyside.*Sin;                                                                                                                     % Heavyside function kills off Sin for negative Sin (energy and momentum transferring from waves to wind) and for negative/zero energy

           % Set Sin = 0 on land boundaries to avoid newdelt --> 0.
           Sin(D<=0) = 0;
           
           p = model.p;
           for tj = 1:p
               short(:,:,:,tj) = sum(E.*cth2(:,:,:,rem((1:p)-tj+p,p)+1),4)*dth;                                                                      % energy in each angular bin get the mean square slope (eqn. 16, Donelan 2012)
           end
           short = (cumsum((wn.^3.*short.*dwn),3)-wn.^3.*short.*dwn);                                                                                % sqrt of mean square slope (eqn. 16, Donelan 2012)
          
% -- Sdt ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           % Calculate turbulent dissipation: Sdt = 4 nu_t k^2.
           Sdt(:,:,:,:) = repmat(Ustar,[1 1 model.o model.p]);
           Sdt = Sdt_fac*sqrt(rhorat).*Sdt.*wn;                                                                                                 % eqn 20, Donelan 2012
%-- Sbf -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           Sbf = zeros(model.m,model.n,model.o,model.p);
           Sbf(D>0) = Sbf_fac*wn(D>0)./sinh(2*wn(D>0).*D(D>0));                                                                                 % eqn. 22, Donelan 2012
% -- Sds ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           Sds = zeros(size(Sin));
           Sds(D>0)=abs(ann(D>0)*2*pi.*f(D>0).*(1+mss_fac*short(D>0)).^2.*(wn(D>0).^4.*E(D>0)).^(nnn(D>0)));                                    % LH p 712. vol 1 [eqn. 17, Donelan 2012]
          
           % Set Sds = 0 on land boundaries to avoid newdelt --> 0.
           Sds(D<=0) = 0;
          
           % Spread Snl to 2 next longer wavenumbers exponentially decaying as distance from donating wavenumber.
           Snl = zeros(model.m,model.n,model.o,model.p);
           Snl(:,:,1:end-1,:) = bf1*Snl_fac*(Sds(:,:,2:end,:).*E(:,:,2:end,:).*wn(:,:,2:end,:).*dwn(:,:,2:end,:));                              % eqn. 21, Donelan 2012 (first part of sum)
           % a quantity of energy proportional to energy dissipated is passed to longer waves in next two lower wn bins
           Snl(:,:,1:end-2,:) = Snl(:,:,1:end-2,:) + bf2*Snl_fac*(Sds(:,:,3:end,:).*E(:,:,3:end,:).*wn(:,:,3:end,:).*dwn(:,:,3:end,:));         % eqn. 21, Donelan 2012 (second part of sum)
          
           % Renormalize to receiving wavenumber and bandwidth.
           Snl(:,:,:,:) = Snl(:,:,:,:)./(wn(:,:,:,:).*dwn(:,:,:,:));
           % Remove downshifted energy
           Snl(:,:,:,:) = Snl(:,:,:,:) - Snl_fac*(Sds(:,:,:,:).*E(:,:,:,:));                                                                    % eqn. 21, Donelan 2012 (last part of sum)
           Snl(D <= 0) = 0;
% -- Sds_wc ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           % Integrate source functions
           Sds_wc = Sds;                                                                                                         % Keep whitecapping (wc) only dissipation for calculating snl growth.
           Sds(D>0) = coth(0.2*wn(D>0).*D(D>0)).*Sds(D>0) + Sdt(D>0) + Sbf(D>0) + 4*liquid.nu_liquid *wn(D>0).^2;                % Add viscous, turbulent and plunging dissipation after calculation of Snl
           % aa = input - dissipation
           aa = Sin(:,:,ol,:) - Sds(:,:,ol,:);                                                                                   % ol = long wavelength waves
           aa = aa(D(:,:,ol,:)>0);                                                                                                    
           aa = max(abs(aa(:)));
           aaa = explim;
          
           if isnan(aa)                                                                                                          % if source = dissipation then denominator for new possible time step is 5e-5
               aa = aaa/model.maxdelt;
           end
           newdelt = max([aaa/aa model.mindelt]);                                                                                % newdelt = delt to give max of 50% growth or decay
           newdelt = min([newdelt model.maxdelt (model.time_step-sumt) delt]);                                                   % min[max(0.1/(Sin-Sds) 0.0001) 2000 TotalModelTime CourantGrid]
           fprintf('newdelt: %.2f\n',newdelt);
          
           % add to model time
           sumt = sumt + newdelt;
           modt = modt + newdelt;
          
           fac_exp = (Sin(:,:,ol,:) - Sds(:,:,ol,:));
          
% -- Wave Energy ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           E1 = zeros(size(E));                                                                                               % E^(n+1) in forward Euler differencing
           E1(:,:,ol,:) = E(:,:,ol,:).*exp(newdelt*fac_exp);                                                                  % Long waves.
          
          
           E1(:,:,ol,:) = E1(:,:,ol,:) + newdelt*Snl(:,:,ol,:);                                                               % eqn. B3, Donelan 2012 (time discretizaton solition for variance spectrum at next time step)

           cath = ones(size(Sds));                                                                                            % horizontal-to-vertical orbital velocity enhancement which leads to more rapid dissipation in shoaling waves relative to deep water spilling breakers
           cath(D>0) = coth(0.2*wn(D>0).*D(D>0));                                                                             % limits the breaker height to depth of shoaling wave ratio [eqn. 17, Donelan 2012 (but A2 = 42 not 0.2?)  
           
           E2(:,:,:,:) = zeros(model.m,model.n,model.o,model.p);
           fij = find((Sin(:,:,:,:) - 4*liquid.nu_liquid*wn(:,:,:,:).^2 - Sdt(:,:,:,:) - Sbf(:,:,:,:)) > 0);                  % finds terms where input > dissipation
           E2(fij) = wn(fij).^(-4).*((Sin(fij)-4*liquid.nu_liquid*wn(fij).^2 - Sdt(fij) - Sbf(fij))./(cath(fij) ...
               .*ann(fij)*2*pi.*f(fij).*(1+mss_fac*short(fij)).^2)).^(nnninv(fij));                                           % LH p.712, vol 1
           E2(D<=0) = 0; E1(D<=0) = 0; 
           E2(E2<0) = 0; E1(E1<0) = 0;
          
           E1(:,:,os,:) = E2(:,:,os,:);                                                                                       % E1 = E2 at small wavelengths (os = short frequency)                                                                                      
          
           E1 = real(E1); E1(E1 < 0)=0;E1(D <= 0)=0;
           E(D <= 0) = 0; E(:,:,os,:) = E1(:,:,os,:);
           E(:,:,ol,:) = (E(:,:,ol,:)+E1(:,:,ol,:))/2;                                                                       % E = (E+E1)/2 (time-splitting for stable integration) [eqn. B5, Donelan 2012]
          
% -- Compute advection term -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           advect = zeros(model.m,model.n,model.o,model.p);
           Ccg = Cg + Uer.*cth + Uei.*sth;
          
           % Upwave advection with account taken of varying delx and dely
           advect(:,:,ol,cp) = advect(:,:,ol,cp)+(Ccg(:,:,ol,cp).*cth(:,:,ol,cp).*E(:,:,ol,cp).*dely(:,:,ol,cp)...                % advection term in eqn. B6, Donelan 2012
               - Ccg(xp,:,ol,cp).*cth(xp,:,ol,cp).*E(xp,:,ol,cp).*dely(xp,:,ol,cp))...
               ./((dely(xp,:,ol,cp) + dely(:,:,ol,cp)).*(delx(xp,:,ol,cp) + delx(:,:,ol,cp)))*4;
          
           advect(:,:,ol,cm) = advect(:,:,ol,cm)+(Ccg(xm,:,ol,cm).*cth(xm,:,ol,cm).*E(xm,:,ol,cm).*dely(xm,:,ol,cm)...            % advection term in eqn. B6, Donelan 2012
               - Ccg(:,:,ol,cm).*cth(:,:,ol,cm).*E(:,:,ol,cm).*dely(:,:,ol,cm))...
               ./((dely(:,:,ol,cm) + dely(xm,:,ol,cm)).*(delx(:,:,ol,cm) + delx(xm,:,ol,cm)))*4;
          
           advect(:,:,ol,sp) = advect(:,:,ol,sp)+(Ccg(:,:,ol,sp).*sth(:,:,ol,sp).*E(:,:,ol,sp).*delx(:,:,ol,sp)...                % advection term in eqn. B6, Donelan 2012
               - Ccg(:,yp,ol,sp).*sth(:,yp,ol,sp).*E(:,yp,ol,sp).*delx(:,yp,ol,sp))...
               ./((dely(:,yp,ol,sp) + dely(:,:,ol,sp)).*(delx(:,yp,ol,sp) + delx(:,:,ol,sp)))*4;
          
           advect(:,:,ol,sm) = advect(:,:,ol,sm)+(Ccg(:,ym,ol,sm).*sth(:,ym,ol,sm).*E(:,ym,ol,sm).*delx(:,ym,ol,sm)...            % advection term in eqn. B6, Donelan 2012
               - Ccg(:,:,ol,sm).*sth(:,:,ol,sm).*E(:,:,ol,sm).*delx(:,:,ol,sm))...
               ./((dely(:,:,ol,sm) + dely(:,ym,ol,sm)).*(delx(:,:,ol,sm) + delx(:,ym,ol,sm)))*4;
          
% -- Full energy from source, sink, and advection -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           E1(:,:,ol,:) = E1(:,:,ol,:) - newdelt*advect(:,:,ol,:);                                                              % eqn. B6, Donelan 2012
          
           % clean up
           E1(D <= 0) = 0;                                                                                                      % energy on land is set to zero
           E1=real(E1);                                                                                                      
           E1(E1 < 0)=0;
          
% -- Compute refraction including wave-current interaction ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
          
           Crot(:,:,ol,:) = (c(xm,:,ol,:)-c(xp,:,ol,:)).*sth(:,:,ol,:)./(delx(xm,:,ol,:)+delx(xp,:,ol,:))...                                   % eqn. A8, Donelan 2012
               -(c(:,ym,ol,:)-c(:,yp,ol,:)).*cth(:,:,ol,:)./(dely(:,ym,ol,:)+dely(:,yp,ol,:));
          
           % Determine if rotation is clockwise or counter clockwise
           Crotcw = zeros(size(Crot));Crotccw = zeros(size(Crot));                  
           Crotccw(Crot>0) = Crot(Crot>0); Crotcw(Crot<0) = Crot(Crot<0);
           Crotccw(Crotccw > dth/newdelt) = dth/newdelt;
           Crotcw(Crotcw < -1*dth/newdelt) = -1*dth/newdelt;

           E1(:,:,ol,cw) = E1(:,:,ol,cw) - newdelt*Crotcw(:,:,ol,:).*E(:,:,ol,:)/dth;                                                          % eqn. A1, Donelan 2012 (clockwise rotation)
           E1(:,:,ol,:) = E1(:,:,ol,:) + newdelt*Crotcw(:,:,ol,:).*E(:,:,ol,:)/dth;                                                            % not sure
           E1(:,:,ol,ccw) = E1(:,:,ol,ccw) + newdelt*Crotccw(:,:,ol,:).*E(:,:,ol,:)/dth;                                                       % eqn. A1, Donelan 2012 (counterclockwise rotation)
           E1(:,:,ol,:) = E1(:,:,ol,:) - newdelt*Crotccw(:,:,ol,:).*E(:,:,ol,:)/dth;                                                           % not sure
           E1(E1 < 0) = 0; E1(D <= 0) = 0;
           E = real(E1);E(isnan(E)) = 0;

           % Make all land points = 0 and upstream values equal to coarse grid
           E(D <= 0) = 0;
                    
           % Compute form stress spectrum.
           tauE = sum(wn.*E.*Sin.*cth./c,4)*dthd*dr;                                                                                            % eastward wind stress (eqn. 5, Donelan 2012)
           tauN = sum(wn.*E.*Sin.*sth./c,4)*dthd*dr;                                                                                            % northward wind stress (eqn. 5, Donelan 2012)
           % Add wind speed dependent tail of slope, "mtail" pinned to the highest wavenumber, "wnh".
           mtail = 0.000112*U_10.*U_10 - 0.01451.*U_10 - 1.0186;
           wnh = squeeze(wn(:,:,model.o,1));                                                                                                    % largest wavenumber
          
           tauE = liquid.rho_liquid*(sum((planet.gravity+liquid.surface_tension.*squeeze(wn(:,:,:,1)).^2./liquid.rho_liquid).*squeeze(dwn(:,:,:,1)).*tauE,3) + ...
               (planet.gravity+liquid.surface_tension.*squeeze(wn(:,:,model.o,1)).^2./liquid.rho_liquid).*squeeze(tauE(:,:,model.o)).*wnh.^(-mtail).*(kutoff.^(mtail+1)-wnh.^(mtail+1))./(mtail+1));    % eqn. 5, Donelan 2012
           tauN = liquid.rho_liquid*(sum((planet.gravity+liquid.surface_tension.*squeeze(wn(:,:,:,1)).^2./liquid.rho_liquid).*squeeze(dwn(:,:,:,1)).*tauN,3) + ...
               (planet.gravity+liquid.surface_tension.*squeeze(wn(:,:,model.o,1)).^2./liquid.rho_liquid).*squeeze(tauN(:,:,model.o)).*wnh.^(-mtail).*(kutoff.^(mtail+1)-wnh.^(mtail+1))./(mtail+1));    % eqn. 5, Donelan 2012
          
           Cd = abs(tauE + i*tauN)./rhoa./(U_z.^2);                                                                                             % shear stress and law of the wall
           Cdf = Cd;
           Ustar_smooth = smooth_nu((1-wfac)*U_z(:),z,planet.nua);
           Ustar_smooth = reshape(Ustar_smooth,model.m,model.n);                                                                                % Surface current (friction velocity for hydraulically smooth interface)
          
           Ustar_smooth = (Ustar_smooth./U_z).^2;
           Cds =  Ustar_smooth;
          
           Ustar_smooth = U_z.^2.*(0.3333*Ustar_smooth + 0.6667*(Ustar_smooth.^2)./(Ustar_smooth + Cd));
           tauE = tauE + rhoa*Ustar_smooth.*cos(squeeze(windir(:,:,1,1)));                                                                      % wind stress + wind momentum in Eastward direction
           tauN = tauN + rhoa*Ustar_smooth.*sin(squeeze(windir(:,:,1,1)));                                                                      % wind stress + wind momentum in Northward direction

           Cd = abs(tauE + i*tauN)./rhoa./(U_z.^2);                                                                                             % form drag coefficient (eqn. 8, Donelan 2012)
          
           
           if Etc.showplots && rem(tplot,10) == 0                                                                                               % plot every 10th time step if showplots = 1
              
               close all;

               % diagnostic plot
               figure(16);clf;semilogx(squeeze(wn(2,model.lati,:,model.p/2)),squeeze(sum(wn(2:4:model.m,model.lati,:,:).*E(2:4:model.m,model.lati,:,:),4)*dthd*dr)','*-');grid on
               title(['Omni directional wavenumber spectra along latitude ',num2str(model.lati),'. Time = ',num2str(modt/3600),' Hours'])      
               % Drag coefficient vs fetch
               figure(202);clf;subplot(321);plot(Cd(:,model.lati).*(U_z(:,model.lati)./U_10(:,model.lati)).^2,'.-');
               title('Drag coefficient')
               grid on
               figure(202);hold on;subplot(322);semilogx(squeeze(wn(model.long,model.lati,:,10)),squeeze(sum(wn(model.long,model.lati,:,:).^2.*Sin(model.long,model.lati,:,:).*E(model.long,model.lati,:,:),4))*dthd*dr,'*-');
               title('k*Sin');grid on
               figure(202);hold on;subplot(323);semilogx(squeeze(wn(model.long,model.lati,:,10)),squeeze(sum(wn(model.long,model.lati,:,:).^2.*Sds(model.long,model.lati,:,:).*E(model.long,model.lati,:,:),4))*dthd*dr,'*-');
               figure(202);hold on;subplot(323);semilogx(squeeze(wn(model.long,model.lati,:,10)),squeeze(sum(wn(model.long,model.lati,:,:).^2.*Sdt(model.long,model.lati,:,:).*E(model.long,model.lati,:,:),4))*dthd*dr,'*-g');
               figure(202);hold on;subplot(323);semilogx(squeeze(wn(model.long,model.lati,:,10)),squeeze(sum(wn(model.long,model.lati,:,:).^2.*Sbf(model.long,model.lati,:,:).*E(model.long,model.lati,:,:),4))*dthd*dr,'*-r');
               title('k*Sds(b), k*Sdt(g), k*Sbf(r)');grid on
               figure(202);hold on;subplot(324);semilogx(squeeze(wn(model.long,model.lati,:,10)),squeeze(sum(wn(model.long,model.lati,:,:).^2.*Snl(model.long,model.lati,:,:),4))*dthd*dr,'*-');
               title('k*Snl');grid on
               figure(202);hold on;subplot(325);semilogx(squeeze(wn(model.long,model.lati,:,10)),squeeze(sum(wn(model.long,model.lati,:,:).^2.*E(model.long,model.lati,:,:),4))*dthd*dr,'*-');
               title('k*spectrum');grid on
           end

% -- Sig wave height -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
          % integrate spectrum to find significant height  
                ht = sum(dwn.*wn.*E,4)*dthd*dr;                            % integral = sum(wavespectrum*wn*del_wn) -->  spectral moment
                ht = sum(ht,3);                                            % sum the prev sum over all frequencies to get the zeroth order moment (aka variance of sea surface (1/2a^2))
                ht = 4*sqrt(abs(ht));                                      % signifigant wave height (from zeroth order moment of surface)
                [sigH(t),~] = max(max(ht));                            % return the largest signifigant wave height along the grid
                htgrid = ht;                                          % return signifigant wave height at each spatial point (m,n) on the grid
              

            if Etc.showplots && rem(tplot,10) == 0    
               % Plot Signifigant Height
               [xplot,yplot] = meshgrid(1:model.m,1:model.n);
               figure;
               surf(xplot',yplot',ht,'EdgeColor','k')
               myc = colorbar;
               myc.Label.String = 'Sig H [m]';
               title(['Sig Wave Height for u = ',num2str(wind.speed),' m/s'])
               frame = getframe(gcf);
               im{idx} = frame2im(frame);
               idx = idx + 1;
            end
% -- mean slope ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           % integrate spectrum to find mean slope                                                                                                                         
            ms = sum(dwn.*wn.^3.*E,4)*dthd*dr;                                                                                                              % slope = angle water surface makes with flat surface
            ms = sum(ms,3);
            ms = sqrt(ms);                                                                                                                                  % standard deviation of water surface = sqrt(variance of water surface)

           if Etc.showplots && rem(tplot,10) == 0 
               figure(202);hold on;subplot(326);plot(1:model.m,ht(:,model.lati),'.-',1:model.m,ms(:,model.lati),'--r');
               title('Sig. Ht. & mean slope');
               grid on
               pause(3)
               figure(20);clf;pcolor(ht');shading interp;colorbar;grid
           end

           % Integrate spectrum over wavenumber to plot directional plot of 10th wavelength.
           Wavel = sum(squeeze(wn(:,:,10,:)).*squeeze(E(:,:,10,:).*dwn(:,:,10,:)),3)./sum(squeeze(E(:,:,10,:).*dwn(:,:,10,:)),3);                           % not sure
           Wavel = 2*pi./Wavel;
           Wavel = Wavel.^(0.25); 
           % Direction of one wavelength.
           KD = squeeze(sum(dwn(:,:,10,:).*E(:,:,10,:),3));                                                                                                 % direction vector of wave at 10th frequency bin (near center for o = 25)
           Dir_sin = sum(KD.*squeeze(sin(waveang(:,:,10,:))),3)./sum(KD,3);                                                                                 % x-component of wave at 10th frequency bin
           Dir_cos = sum(KD.*squeeze(cos(waveang(:,:,10,:))),3)./sum(KD,3);                                                                                 % y-component of wave at 10th frequency bin
           mdir = atan2(Dir_sin,Dir_cos);                                                                                                                   % average wave direction [-pi to pi]

           if Etc.showplots && rem(tplot,10) == 0
               figure(20);
               hold on;
               quiver(Wavel'.*cos(mdir'),Wavel'.*sin(mdir'),0.8,'c');
               title(['Sig.Ht. and \lambda^{1/4}. Time = ',num2str(modt/3600),' Hours'])
               pause(2)
               figure(36);
               clf;
               contour(ht',20);
               colorbar;
               grid on
               title(['Sig.Ht. Time = ',num2str(modt/3600),' Hours'])
           end
          
           cgmax = max(max(max(max((E>1e-320).*Cg))));                                                                                                      % Adjust cgmax for active energy components.
          

       end % end while loop for sub-time steps
      
       fraction_time_completed = t/model.num_time_steps;
       fprintf('u = %.2f: fraction time completed: %.2f\n',UU,fraction_time_completed);
    
      E_each{t} = E;                                                                                                                                    % return energy spectrum for wind speed at each time step t

      if Etc.savedata
        save([TitanResults,'/New_u_',int2str(UU),'_t_',int2str(t),'_',Etc.name],'E','ht','freqs','oa','Cd','Cdf','Cds','Sds','Sds_wc','Sin','Snl','Sdt','Sbf','ms')
      end
       
       
       if t > 100 && sigH(t-1)/sigH(t) < model.tolH                                                                                                 % will break out of wind speed loop if waves haven't changed by more than the tolerance level tolH 
           disp('Waves have reached 99% of maturity.')
           break
       elseif t > 1
           fprintf('t_n-1/t_n: %.6f\n',sigH(t-1)/sigH(t));
       end

   end % end of loop for Tsteps
  
   [~,minind] = min(abs(wind.speed-UU));
   disp(['Finished Wind Speed ' num2str(UU) ' m/s (' num2str(minind) ' Of ' num2str(numel(wind.speed)) ' )']);
end
