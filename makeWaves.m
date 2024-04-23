function [sigH,htgrid,E_each,ms] = makeWaves(planet,model,wind,uniflow,Etc)
%% ==========================================================================================================================================================================================================================================================================
%% ==========================================================================================================================================================================================================================================================================
% MAKEWAVES calculates E(x,y,k,theta) for wave field using an energy balance between wind input and multiple dissipation terms including turbulent dissipation (Sdt), bottom friction (Sbf), wave breaking (Sds), and spilling breakers (Ssb) as well as a non-linear
% interaction term that shifts energy conservatively within the wave spectrum
%
% The equation to solve is:
%   E_{n+1} = E_{n} + del*(-Cg*cos(theta)*dE/dx - Cg*sin(theta)*dE/dy] + Sin + Snl - Sds
%
%
% Dimensions of (x,y,k,theta) are (m,n,o,p)
%
%
%   Arguments:
%       planet
%           rho_liquid      : liquid density [kg/m3]
%           nu_liquid       : liquid kinematic viscocity [m2/s]
%           nua             : atmospheric gas viscocity [m2/s]
%           gravity         : gravitational acceleration [m/s2]
%           surface_temp    : surface temperature [K]
%           surface_press   : surface pressure [Pa]
%           surface_tension : surface tension of liquid [N/m]
%           name            : planet name ['string'] e.g. 'Titan'
%       model
%           m               : number of grid cells in x-dimension
%           n               : number of grid cells in y-dimension
%           o               : number of frequency bins
%           p               : number of angular (theta) bins. Must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
%           long            : longitude location of grid to plot 
%           lat             : latitude location of grid to plot
%           gridX           : size of grid cell in x-dimension [m]
%           gridY           : size of grid cell in y-dimension [m]
%           mindelt         : minimum time step to minimize oscillations of delt to very small values [s]
%           maxdelt         : maximum time step to prevent large values as wind ramps up at startup [s]
%           tolH            : minimum change in SigH to stop model at (otherwise runs to full Tsteps) 
%           cutoff_freq     : cutoff frequency bin seperating the diagnostic from advecting wavenumbers
%           min_freq        : minimum frequency to model [Hz]
%           max_freq        : maximum frequency to model [Hz]
%           bathy_map       : m x n array of depth [m] (+ values = subsurface, - values = subaerial)
%           time_step       : maximum size of time step [s]
%           num_time_steps  : length of model run in terms of # of time steps
%           tune_A1         : tuning parameter for input term (Donelan et al. 2012, eqn. 12)
%           tune_mss_fac    : A3 tuning parameter in spilling breaker term (Donelan et al. 2012, eqn. 16)
%           tune_Sdt_fac    : fraction of Sdt term going into spectrum (A4 in Donelan et al. 2012, eqn. 20)
%           tune_Sbf_fac    : fraction of Sbf term going into spectrum (Gf in Donelan et al. 2012, eqn. 22)
%           tune_cotharg    : tuneable constant in argument in coth() term for dissipation term (is 1 but later improved fit using 0.2 in Donelan et al. 2012, eqn. 17)
%           tune_n          : exponent n term in dissipation term (Donelan et al. 2012, eqn. 15)
%           z_data          : elevation where wind measurements are taken (e.g. u10 would correspond to z_data = 10 m) [m]
%       wind
%           dir             : direction of incoming near-surface wind (CCW from East) [radians]
%           speed           : magnitude of incoming near-surface wind [m/s]
%       uniflow
%           East            : Eastward unidirectional current [m/s]
%           North           : Northward unidirectional current [m/s]
%       Etc
%           showplots       : 0 = no plots made, 1 = plots every tenth loops (will slow down model run)
%           savedata        : 1  = save data for time steps (will slow down model run), 0 = skip saving
%           showlog         : 1 = print progress to command line (will slow down model run), 0 = no progress printed to command line
%
%
%   Returns:
%       sigH                : significant wave height at specified (Model.lat, Model.lon) coordinates [m]
%       htgrid              : significant wave height for each grid cell [m]
%       E_each              : wave energy spectrum (x,y) in space and in (frequency,direction) space
%       ms                  : mean slope of liquid surface
%     
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% External function requirements: none
% Compatible with MATLAB 2024a
%
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Approx time to run: 7 minutes (for 1 wind speed, 10 time steps for uniform depth)
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Cite: 
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% Authors: Mark Donelan, Alex Hayes, Charlie Detelich, Una Schneck
tic
%% ==========================================================================================================================================================================================================================================================================
%% ==========================================================================================================================================================================================================================================================================

% -- save input directory to log --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TitanResults = make_log(planet,model,wind,uniflow,Etc);

% -- prepare .mat files to be saved during loops -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
file = 1;

% -- MODEL SET UP ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% global constants
kappa = 0.4;                                                               % Von-Karman constant
i = sqrt(-1);                                                              % imaginary number i
kgmolwt = 0.028;                                                           % gram molecular weight [Kgm/mol]
RRR = 8.314;                                                               % Universal gas constant [J/K/mol]
% frequency and direction bins
dr = pi/180;                                                               % conversion from degrees to radians
dthd = 360/(model.Dirdim);                                                      % step size for angular direction [degrees]
dth = dthd*dr;                                                             % 2pi/p (small angle for integration [radians])
% Set up geographic deltas
dely = model.gridY*ones(model.LonDim,model.LatDim,model.Fdim,model.Dirdim);
delx = model.gridX*ones(model.LonDim,model.LatDim,model.Fdim,model.Dirdim);
% Wind:
gust = 0;                                                                  % Gust factor applied to wind at each time step
zref = 20;                                                                 % height of reference wind speed, normally at 20m [m]
wfac = 0.035;                                                              % wind drift fraction of Uz ffor U10m
% Ideal gas law: PV = nRT
% Densities:
rhoa = planet.surface_press*kgmolwt/(RRR*planet.surface_temp);             % air density [kg/m3]
rhorat=rhoa/planet.rho_liquid;                                             % air-water density ratio.
% Wavenumber limits:
kutoff = 1000;                                                             % wavenumber cut-off due to Kelvin-Helmholtz instabilities (Donelan 2012, pg. 3)
kcg = sqrt(planet.gravity*planet.rho_liquid/planet.surface_tension);       % wavenumber of slowest waves defined by the capillary-gravity waves, from Airy dispersion: omega^2 =gktanh{k)
% modified wavenumbers
kcga = 1.15*kcg;                                                           % a min is shifted above kcg.
% create frequency limits for spectrum
dlnf=(log(model.max_freq)-log(model.min_freq))/(model.Fdim-1);                % frequency step size for log normal distribution
f = exp(log(model.min_freq)+(0:model.Fdim-1)*dlnf);                           % frequencies for spectrum
dom = 2*pi*dlnf.*f;                                                        % discrete angular frequency (w = 2pi*f)
freqs = f;                                                                 % save a copy of frequencies

if Etc.showplots
    figure;
    plot(1:numel(freqs),freqs,'-','LineWidth',3)
    hold on
    xline(model.cutoff_freq,'-','Cutoff Frequency')
    xlabel('frequency bin')
    ylabel('frequency [hz]')
    set(gca, 'YScale', 'log')
end

% Frequency bins:
ol = 1:model.cutoff_freq;                                                  % bins for long frequencies
os = model.cutoff_freq+1:model.Fdim;                                          % bins for short frequencies
% Compute diffusion values in 2 freqs and 2 directions:
%   bf1 + bf2 = 1, and bt1 + bt2 = 1.
bfac = 16;                                                                 % has also been set to 15.73 by Donelan in past
bf1 = exp(-bfac*dlnf*dlnf);                                                % part of non-linear source term for downshifting and spilling breakers [eqn. 21 Donelan+2012]
bf2 = exp(-bfac*4*dlnf*dlnf);                                              % part of non-linear source term for downshifting and spilling breakers [eqn. 21 Donelan+2012]
bf1a = bf1/(bf1 + bf2);                                                    % normalization [eqn. 21, Donelan+2012]
bf2 = bf2/(bf1 + bf2);                                                     % normalization [eqn. 21, Donelan+2012]
bf1 = bf1a;
% Compute Snl_fac
A = exp(dlnf);
B = A.*A;
fac = bf1.*(1-1/B) + bf2.*(1-1./B.^2);                                     % eqn. 21, Donelan+2012
Snl_fac = ((1.0942/A)^1.9757)/fac;                                         % eqn. 21, Donelan+2012
% waveangle [radians]
waveang = ((0:model.Dirdim-1)-model.Dirdim/2+0.5)*dth;                               % wave angle [Radians]
th = ((0:model.Dirdim-1)-model.Dirdim/2+0.5)*dth;                                    % angle of wave propagation (phi in Donelan+2012)  between -pi to +pi
cth=cos(th);                                                               % cosine of angle in wave propagation direction 
sth=sin(th);                                                               % sine of angle in wave propagation direction
% compute cos^2 for calculation of mss vs direction.
tth = (0:model.Dirdim-1)*dth;                                                   % angular difference between short waves and longer waves
cth2 = cos(tth).^2;                                                        % cosine square of angular difference between short waves and longer waves
% indices for refraction rotation:
cw = ([model.Dirdim 1:model.Dirdim-1]);                                              % clockwise rotating indices
ccw = ([2:model.Dirdim 1]);                                                     % counterclockwise rotating indices
% upwave indices for advective term:
xp = [1 1:model.LonDim-1];
yp = [1 1:model.LatDim-1];
xm = [2:model.LonDim model.LonDim];
ym = [2:model.LatDim model.LatDim];
% Index being used for advection:
cp = find(cth > 0);
cm = find(cth < 0);
sp = find(sth > 0);
sm = find(sth < 0);
% reshape matrices to 4D by repeating over all other dimensions
cth = repmat(cth,[model.Fdim 1 model.LonDim model.LatDim]);
cth=shiftdim(cth,2);
sth = repmat(sth,[model.Fdim 1 model.LonDim model.LatDim]);
sth=shiftdim(sth,2);
cth2 = repmat(cth2,[model.Fdim 1 model.LonDim model.LatDim]);
cth2=shiftdim(cth2,2);
%% -- Lake Geometry ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
D = model.bathy_map;                                                       % depth of liquid [m]
D(D<=0) = 0;                                                               % limit land elevations to 0 to avoid dD/dx, dD/dy errors in refraction calculation    
D(1, :) = 0;D(end, :) = 0;D(:, 1) = 0;D(:, end) = 0;                       % set depth array boundary to 0 (absorptive boundary condition)
% plot the bathymetry
[xplot,yplot] = meshgrid(1:model.LonDim,1:model.LatDim);
if Etc.showplots
    figure;
    h1 = surf(xplot,yplot,D','EdgeColor','k','FaceColor','interp','FaceAlpha',0.5);
    myc = colorbar;
    myc.Label.String = 'Liquid Depth [m]';
    title('Lake Model Bathymetry')
    view(2)
    hold on
    h2 = quiver(model.long,model.lat,cos(wind.dir),sin(wind.dir),'r','LineWidth',2,'MaxHeadSize', 1);
    legend(h2,'Wind Direction', 'Location', 'northwest','interpreter','latex');
    hold off
end
%% -- wavemumber and Power of Sds ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% initalize 4D arrays to be filled 
wn = ones(model.LonDim,model.LatDim,model.Fdim);
nnn = ones(model.LonDim,model.LatDim,model.Fdim);
ann = ones(model.LonDim,model.LatDim,model.Fdim);

for ii = 1:model.LatDim
   for jj = 1:model.LonDim
       if D(ii,jj) > 0
           wn(jj,ii,:) = wavekgt(f,D(ii,jj),planet.gravity,planet.surface_tension,planet.rho_liquid,1e-2);                   % wave number (using linear wave dispersion)
           nnn(jj,ii,:) = model.tune_n;                                                                                      % Power n of Sds on the degree of saturation [eqn. 15, Donelan+2012, eqn. 5, Donelan+2001] 1.2 + 1.3*(abs(2 - (1+3*(wn(jm,jn,:)./kcgn).^2)./(1+(wn(jm,jn,:)./kcgn).^2)).^2.0)
           ann(jj,ii,:) = 0.04 + 41.96*(abs(2 - (1+3*(wn(jj,ii,:)./kcga).^2)./(1+(wn(jj,ii,:)./kcga).^2)).^4.0);             % Power of Sds
       else
           wn(jj,ii,:) = NaN(model.Fdim,1);
           nnn(jj,ii,:) = NaN(model.Fdim,1);
           ann(jj,ii,:) = NaN(model.Fdim,1);
       end
   end
end

% reshape the matrix by repeating over the angle dimension p (dimensions of [n m o p])
wn = repmat(wn,[1 1 1 model.Dirdim]);
nnn = repmat(nnn,[1 1 1 model.Dirdim]);
ann = repmat(ann,[1 1 1 model.Dirdim]);

if Etc.showplots

    surf_extrema(wn,'wn',model,'tallest')
    surf_extrema(wn,'wn',model,'smallest')
    plot_freq_depend(wn,'wn',D,freqs,model) 


    surf_extrema(nnn,'n',model,'tallest')
    surf_extrema(nnn,'n',model,'smallest')
    plot_freq_depend(nnn,'n',D,freqs,model)

    surf_extrema(ann,'an',model,'tallest')
    surf_extrema(ann,'an',model,'smallest')
    plot_freq_depend(ann,'an',D,freqs,model)

end

% conversions for nnn?
nnn = (2.53/2.5).*nnn;
nnninv = 1./nnn;
% reshape the matrix by repeating over the angle dimension p (dimensions of [n m o p])
D = repmat(D',[1 1 model.Fdim model.Dirdim]);                                                                                         % 2D matrix of depths repeated in the frequency and direction dimension (D(x,y,i,j) = D(x,y) for all i and j)
f = repmat(f',[1 model.Dirdim model.LonDim model.LatDim]);                                                                                  % 1D matrix of frequencies repeated in x, y, and direction dimension (f(i,j,A,k) = f(A) for all i,j,k)
f = shiftdim(f,2);
dom = repmat(dom',[1 model.Dirdim model.LonDim model.LatDim]);                                                                              % 1D matrix of discrete angular frequencies repeated in x, y, and direction dimension (dom(i,j,A,k) = dom(A) for all i,j,k)
dom = shiftdim(dom,2);
%% -- ocean currents ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Uer = uniflow.East*D + 0.0;                                                                                                  % Eastward current, m/s
Uei = uniflow.North*D + 0.0;                                                                                                 % Northward current, m/s
%% -- wave speed and group velocity --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
c = (2*pi*f)./wn;                                                                                                            % phase velocity                                                                                                                                                                                            % phase speed
c(D<=0) = 0;                                                                                                                 % set phase speed on land to zero                                                                                        

if Etc.showplots

    surf_extrema(c,'c',model,'tallest')
    surf_extrema(c,'c',model,'smallest')
    plot_freq_depend(c,'c',D,freqs,model)

end

Cg = zeros(size(c));                                                                                                                                                                                                 % initialize group speed
Cg(D>0) = c(D>0)./2.*(1 + 2*wn(D>0).*D(D>0)./sinh(2*wn(D>0).*D(D>0)) + 2*planet.surface_tension.*wn(D>0)./planet.rho_liquid./(planet.gravity./wn(D>0) + planet.surface_tension.*wn(D>0)./planet.rho_liquid));        % Group velocity for all waves (Kinsman)

if Etc.showplots

    surf_extrema(Cg,'Cg',model,'tallest')
    surf_extrema(Cg,'Cg',model,'smallest')
    plot_freq_depend(Cg,'Cg',D,freqs,model)

end

dwn = ones(model.LonDim,model.LatDim,model.Fdim,model.Dirdim);                                                                                                                                                                         % initialize dominant wavenumber (c = dw/dk)
dwn(D>0) = dom(D>0)./abs(Cg(D>0));                                                                                                                                                                                   % remove any values on land for dominant wavenumber (c = dw/dk)

if Etc.showplots

    surf_extrema(dwn,'dwn',model,'tallest')
    surf_extrema(dwn,'dwn',model,'smallest')
    plot_freq_depend(dwn,'dwn',D,freqs,model)

end

% length-scales of interest                                                                                                                                                                                          
l2=abs(c)./(2.*f);                                                                                                                                                                                                  % half the wavelength
lz=abs(c)./(2*pi.*f);                                                                                                                                                                                               % 1/wavenumber 
% Set l2, lz > zref equal to zref   
l2(l2 > zref) = zref;
lz(lz > zref) = zref;

if Etc.showplots

    surf_extrema(l2,'l2',model,'tallest')
    surf_extrema(l2,'l2',model,'smallest')
    plot_freq_depend(l2,'l2',D,freqs,model)

    surf_extrema(lz,'lz',model,'tallest')
    surf_extrema(lz,'lz',model,'smallest')
    plot_freq_depend(lz,'lz',D,freqs,model)

end

E = zeros(model.LonDim,model.LatDim,model.Fdim,model.Dirdim);                                                                                                                                                                          % initializing the energy term in x,y,k,theta with zeros
cgmax = max(max(max(max(Cg))));                                                                                                                                                                                      % fastest group velocity
mindelx = min(squeeze(delx(1,:,1,1)));                                                                                                                                                                               % smallest spatial grid spacing for domain of influence
waveang = repmat(waveang,[model.Fdim 1 model.LonDim model.LatDim]);
waveang=shiftdim(waveang,2);

if Etc.savedata
    m = model.LonDim; n = model.LatDim; o = model.Fdim; p = model.Dirdim; oa = model.cutoff_freq;
    save([TitanResults,'/New_Reference'],'m','n','o','p','freqs','f','wn','dwn','D','Uei','Uer','Cg','c','delx', 'dely','dthd','waveang','dr')
end

%% -- loop through wind speeds --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
idx = 1;                                                                   % frame for making gif
sigH = zeros(1,length(1:model.num_time_steps));                            % initialize sigH for returning (initialize to NaN instead?)
htgrid = cell(1,length(1:model.num_time_steps));                           % initialize htgrid for returning
E_each =  cell(1,length(1:model.num_time_steps));                          % initialize E-spectrum for returning

UU = wind.speed;                                                 

U = UU.*ones(model.LonDim,model.LatDim);                                             % set wind velocity everywhere in x-y plane
windir = wind.dir.*ones(model.LonDim,model.LatDim);                                  % set wind direction everywhere in x-y plane
windir = repmat(windir,[1 1 model.Fdim model.Dirdim]);                             % reshape wind direction matrix by repeating over the frequency and direction arrays in o and p

U_z = U; %+ 0.005;                                                         % wind at modeled height (plus small number to wind speed to avoid division by zero)

% drag coefficient 
Cd = 1.2*ones(size(U));                                                    % drag coefficient for weak/moderate winds 
Cd(U > 11) = 0.49 + 0.065*U(U > 11);                                       % drag coefficient for strong winds (Large&Pond1981, Donelan2004)
Cd = Cd/1000;                                                              % fix units 

modt = 0;                                                                  % model time initialization

% Currents superimposed onto wave-driven flow
Uer = 0*D + 0.0;                                                           % Eastward current [m/s] 
Uei = 0*D + 0.0;                                                           % Northward current [m/s] 

%% -- loop through model time to grow waves --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
disp('================================================================')
disp('starting model time iteration')
disp('================================================================')

ks = eps(1); % start with still surface but can't be exactly zero to define the rough boundary lengthscale

for t = 1:model.num_time_steps                                                                                                                                   % loop through time
  
   sumt = 0;                                                                                                                                                     % initialize total time within time-step t
   tplot = - 1; 
   while ((model.time_step - sumt) > 0)                                                                                                                          % each time step is determined as min[max(0.1/(Sin-Sds) 0.0001) 2000 UserDefinedTime CourantGrid], iterate until the user-defined time is larger than the time passed within time-step t
      
       if gust > 0
            U = U.*(1+gust*randn(size(U)));                                                                                                                      % add random bursts of gusts to wind speed
       end

       tplot = tplot + 1;
       explim = 0.1;                                                                                                                                             % limit for making dynamic time step if the source = dissipation
       delt = 0.7 * mindelx / cgmax;                                                                                                                             % advection-limited Courant condition.
      
% Scale wind to reference heights at 1/wavenumber and wavelength/2 --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------                                                                                                                                           % kinematic stress (u_*^2) at boundary is proportional to square of wind velocity at 10m [u_*^2 = Cd*U^2]  
       Ustar = real(U_z.*sqrt(Cd));                                                                                                                              % frictional velocity of wind at surface               
       U_10 = (1/kappa)*Ustar.*log(10/model.z_data) + U_z;                                                                                                       % log wind profile to scale input wind speed at zdata to speed at 10-m
       ustarw = Ustar .* sqrt(rhorat);                                                                                                                           % scaling on u* for stress calculation later
      
       Ustar = repmat(Ustar,[1 1 model.Fdim model.Dirdim]);
       U = repmat(U_z,[1 1 model.Fdim model.Dirdim]);

       % wind speed scaled to half the wavelength above the surface (function of freq)
       Ul = (1/kappa).*Ustar.*log(l2/model.z_data) + U;                                                                                                           % U(z2)/U(z1) = ln(z2/z0)/ln(z1/z0)
       
       if Etc.showplots% && sum(sum(sum(sum(any(Ul<0))))) > 0 && sumt > 0  
           close all
           plot_freq_depend(Ul,'U_{L/2}',D,freqs,model)
           plot_freq_depend(Ustar,'U_*',D,freqs,model)
           plot_freq_depend(log(l2/model.z_data),'log(l2/z)',D,freqs,model)
           %error('makeWaves: Wind velocity at lambda/2 is negative after model has spun up. Try a different frequency range or resolution.')
           laminar_sublayer_thickness = (Ustar.*ks)./planet.nu_liquid;
           z0_smooth = (planet.nu_liquid)./(9.*Ustar);
           z0_rough = ks./30;
           
           plot_freq_depend(log(l2./z0_smooth),'log(l2/z_0,smooth)',D,freqs,model)
           plot_freq_depend(log(l2./z0_rough),'log(l2/z_0,rough)',D,freqs,model)
           plot_freq_depend(laminar_sublayer_thickness,'laminar sublayer thickness: u*epsilon/nu',D,freqs,model)
            
           pos_range = find(l2>=z0_rough & l2<model.z_data);
           l2_plot = NaN(size(l2));
           l2_plot(pos_range) = l2(pos_range);
           l2_plot = reshape(l2_plot,size(l2));
           test_me_rough = (Ustar./kappa).*log(l2_plot./model.z_data) + U;

           pos_range = find(l2>=z0_smooth & l2<model.z_data);
           l2_plot = NaN(size(l2));
           l2_plot(pos_range) = l2(pos_range);
           l2_plot = reshape(l2_plot,size(l2));
           test_me_smooth = (Ustar./kappa).*log(l2_plot./model.z_data) + U;

           plot_freq_depend(test_me_rough,'(Ustar./kappa).*log(l2_r./model.z_data) + U',D,freqs,model)
           plot_freq_depend(test_me_smooth,'(Ustar./kappa).*log(l2_s./model.z_data) + U',D,freqs,model)


       end
       
       

       ustw = repmat(ustarw,[1 1 model.Fdim model.Dirdim]);
       Ud = - ustw./kappa.*log(lz/model.z_data);                                                                                                                   % drift speed (scaled to 1/wavenumber (function of frequency))
      
       % calculate wind input as a fraction of stress
       relative_angle = windir-waveang;
       Ul_term = Ul.*cos(relative_angle)-c-Ud.*cos(relative_angle)-Uer.*cth-Uei.*sth;
       %_test = Ul.*cos(relative_angle)-c;                                                                                                                         % positive if wind faster than waves, negative if waves faster than wind
    

       Ul_term = abs(Ul_term).*(Ul_term);                                                                                                                          % Sin = A1*[Ul]*((k*wn)/g)*(rhoa/rhow)*E  [eqn. 4, Donelan+2012]
         
       % Adjust input to lower value when waves overrun the wind because as wave speed approaches wind speed, energy extraction becomes less efficient
       Ul_term(Ul_term>0) = model.tune_A1*Ul_term(Ul_term>0);                                                                                                      % wind outruns waves [A1 = 0.11, eqn. 12, Donelan+2012]
       Ul_term(Ul_term<=0) = (model.tune_A1*0.09)*Ul_term(Ul_term<=0);                                                                                             % waves outrun wind [A1 = 0.01, eqn. 12, Donelan+2012]
       Ul_term(cos(windir-waveang)<0) = (model.tune_A1*0.9)*Ul_term(cos(windir-waveang)<0);                                                                        % waves move against wind [A1 = 0.10, eqn. 12, Donelan+2012]

% -- Sin --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       % input to energy spectrum from wind
       Sin = zeros(size(Ul_term));
       Sin(D>0) = rhorat*(Ul_term(D>0).*2*pi.*f(D>0).*wn(D>0)./(planet.gravity+planet.surface_tension.*wn(D>0).*wn(D>0)./planet.rho_liquid));                   % eqn. 4, Donelan+2012
       

       % limits energy going into spectrum once waves outrun wind
       Heavyside = ones(size(E));
       Heavyside(Sin < 0 & E < 1e-320) = 0;

       Sin = Heavyside.*Sin;                                                                                                                                    % Heavyside function kills off Sin for negative Sin (energy and momentum transferring from waves to wind) and for negative/zero energy

       % Set Sin = 0 on land boundaries to avoid newdelt --> 0.
       Sin(D<=0) = 0;
      
       
       for tj = 1:model.Dirdim
           short(:,:,:,tj) = sum(E.*cth2(:,:,:,rem((1:model.Dirdim)-tj+model.Dirdim,model.Dirdim)+1),4)*dth;                                                                   % energy in each angular bin get the mean square slope (eqn. 16, Donelan+2012)
       end
       
       short = (cumsum((wn.^3.*short.*dwn),3)-wn.^3.*short.*dwn);                                                                                               % sqrt of mean square slope [eqn. 16, Donelan+2012]
      
% -- Sdt ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       % Calculate turbulent dissipation: Sdt = 4 nu_t k^2.
       Sdt(:,:,:,:) = Ustar;
       Sdt = model.tune_Sdt_fac*sqrt(rhorat).*Sdt.*wn;                                                                                                         % eqn 20, Donelan+2012

       %-- Sbf -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       Sbf = zeros(model.LonDim,model.LatDim,model.Fdim,model.Dirdim);
       Sbf(D>0) = model.tune_Sbf_fac*wn(D>0)./sinh(2*wn(D>0).*D(D>0));                                                                                         % eqn. 22, Donelan+2012

% -- Sds ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       Sds = zeros(size(Sin));
       Sds(D>0)=abs(ann(D>0)*2*pi.*f(D>0).*(1+model.tune_mss_fac*short(D>0)).^2.*(wn(D>0).^4.*E(D>0)).^(nnn(D>0)));                              % LH p 712. vol 1 [eqn. 17, Donelan+2012]
       % Set Sds = 0 on land boundaries to avoid newdelt --> 0.
       Sds(D<=0) = 0;                                                                                                                            % set dissipation to zero on land

       % Spread Snl to 2 next longer wavenumbers exponentially decaying as distance from donating wavenumber.
       Snl = zeros(model.LonDim,model.LatDim,model.Fdim,model.Dirdim);
       Snl(:,:,1:end-1,:) = bf1*Snl_fac*(Sds(:,:,2:end,:).*E(:,:,2:end,:).*wn(:,:,2:end,:).*dwn(:,:,2:end,:));                                   % eqn. 21, Donelan+2012 (first part of sum)
       % a quantity of energy proportional to energy dissipated is passed to longer waves in next two lower wn bins
       Snl(:,:,1:end-2,:) = Snl(:,:,1:end-2,:) + bf2*Snl_fac*(Sds(:,:,3:end,:).*E(:,:,3:end,:).*wn(:,:,3:end,:).*dwn(:,:,3:end,:));              % eqn. 21, Donelan+2012 (second part of sum)
      
       % Renormalize to receiving wavenumber and bandwidth.
       Snl(:,:,:,:) = Snl(:,:,:,:)./(wn(:,:,:,:).*dwn(:,:,:,:));
       % Remove downshifted energy
       Snl(:,:,:,:) = Snl(:,:,:,:) - Snl_fac*(Sds(:,:,:,:).*E(:,:,:,:));                                                                         % eqn. 21, Donelan+2012 (last part of sum)
       Snl(D <= 0) = 0;                                                                                                                          % set dissipation to zero on land

% -- Sds_wc ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       % Integrate source functions
       Sds_wc = Sds;                                                                                                                             % Keep white-capping (wc) only dissipation for calculating snl growth.
       Sds(D>0) = coth(model.tune_cotharg*wn(D>0).*D(D>0)).*Sds(D>0) + Sdt(D>0) + Sbf(D>0) + 4*planet.nu_liquid *wn(D>0).^2;                     % Add viscous, turbulent and plunging dissipation after calculation of Snl

       % aa = input - dissipation
       aa = Sin(:,:,ol,:) - Sds(:,:,ol,:);                                                                                                       % ol = long wavelength waves (that will advect)
       aa = aa(D(:,:,ol,:)>0);                                                                                                    
       aa = max(abs(aa(:)));
       aaa = explim;
      
       if isnan(aa)                                                                                                                              % if source = dissipation then denominator for new possible time step is 5e-5
           aa = aaa/model.maxdelt;
       end

       newdelt = max([aaa/aa model.mindelt]);                                                                                                    % newdelt = delt to give max of 50% growth or decay
       newdelt = min([newdelt model.maxdelt (model.time_step-sumt) delt]);                                                                       % min[max(0.1/(Sin-Sds) 0.0001) 2000 TotalModelTime CourantGrid]
       %fprintf('newdelt: %.2f\n',newdelt);
      
       % add to model time
       sumt = sumt + newdelt;
       modt = modt + newdelt;
       
       %fprintf('sub-time step sum: %.5f\n',sumt);
       

       fac_exp = (Sin(:,:,ol,:) - Sds(:,:,ol,:));                                                                                                % difference in input and dissipation term to be used in exp() term for energy
      
% -- Wave Energy ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       E1 = zeros(size(E));                                                                                                    % E^(n+1) in forward Euler differencing
       E2 = zeros(size(E));
       E1(:,:,ol,:) = E(:,:,ol,:).*exp(newdelt*fac_exp);                                                                       % Long waves growing with time
      
      
       E1(:,:,ol,:) = E1(:,:,ol,:) + newdelt*Snl(:,:,ol,:);                                                                    % eqn. B3, Donelan 2012 (time discretizaton solution for variance spectrum at next time step)

       cath = ones(size(Sds));                                                                                                 % horizontal-to-vertical orbital velocity enhancement which leads to more rapid dissipation in shoaling waves relative to deep water spilling breakers
       cath(D>0) = coth(0.2*wn(D>0).*D(D>0));                                                                                  % limits the breaker height to depth of shoaling wave ratio [eqn. 17, Donelan 2012 (but A2 = 42 not 0.2?)  (should this 0.2 = Model.cotharg?)
       
       E2(:,:,:,:) = zeros(model.LonDim,model.LatDim,model.Fdim,model.Dirdim);
       fij = find((Sin(:,:,:,:) - 4*planet.nu_liquid*wn(:,:,:,:).^2 - Sdt(:,:,:,:) - Sbf(:,:,:,:)) > 0);                       % finds terms where input > dissipation
       E2(fij) = wn(fij).^(-4).*((Sin(fij)-4*planet.nu_liquid*wn(fij).^2 - Sdt(fij) - Sbf(fij))./(cath(fij) ...
           .*ann(fij)*2*pi.*f(fij).*(1+model.tune_mss_fac*short(fij)).^2)).^(nnninv(fij));                                     % LH p.712, vol 1
       E2(D<=0) = 0; E1(D<=0) = 0; 
       E2(E2<0) = 0; E1(E1<0) = 0;
      
       E1(:,:,os,:) = E2(:,:,os,:);                                                                                            % E1 = E2 at small wavelengths (os = short frequency)                                                                                      
      
       E1 = real(E1); E1(E1 < 0)=0;E1(D <= 0)=0;
       E(D <= 0) = 0; E(:,:,os,:) = E1(:,:,os,:);
       E(:,:,ol,:) = (E(:,:,ol,:)+E1(:,:,ol,:))/2;                                                                             % E = (E+E1)/2 (time-splitting for stable integration) [eqn. B5, Donelan 2012]
      
% -- Compute advection term -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       advect = zeros(model.LonDim,model.LatDim,model.Fdim,model.Dirdim);                                                                        % initialize the 4D advection array
       Ccg = Cg + Uer.*cth + Uei.*sth;                                                                                         % group velocity including unidirectional current and wave velocity
      
       % Upwave advection with account taken of varying delx and dely
       % Note: only long frequencies are advected, short frequencies are assumed to be in equilibrium with
       % the wave and do not need to be advected (model.cutoff_freq defines the cutoff between short and long frequency 
       % bins and can be varied by the user). If not enough frequencies are being advected because of a 
       % poorly chosen model.cutoff_freq, then there will be no fetch dependence across the grid since no energy is coming '
       % into a grid from its neighbor
       advect(:,:,ol,cp) = advect(:,:,ol,cp)+(Ccg(:,:,ol,cp).*cth(:,:,ol,cp).*E(:,:,ol,cp).*dely(:,:,ol,cp)...                % advection term in eqn. B6, Donelan+2012 
           - Ccg(xp,:,ol,cp).*cth(xp,:,ol,cp).*E(xp,:,ol,cp).*dely(xp,:,ol,cp))...
           ./((dely(xp,:,ol,cp) + dely(:,:,ol,cp)).*(delx(xp,:,ol,cp) + delx(:,:,ol,cp)))*4;
      
       advect(:,:,ol,cm) = advect(:,:,ol,cm)+(Ccg(xm,:,ol,cm).*cth(xm,:,ol,cm).*E(xm,:,ol,cm).*dely(xm,:,ol,cm)...            % advection term in eqn. B6, Donelan+2012
           - Ccg(:,:,ol,cm).*cth(:,:,ol,cm).*E(:,:,ol,cm).*dely(:,:,ol,cm))...
           ./((dely(:,:,ol,cm) + dely(xm,:,ol,cm)).*(delx(:,:,ol,cm) + delx(xm,:,ol,cm)))*4;
      
       advect(:,:,ol,sp) = advect(:,:,ol,sp)+(Ccg(:,:,ol,sp).*sth(:,:,ol,sp).*E(:,:,ol,sp).*delx(:,:,ol,sp)...                % advection term in eqn. B6, Donelan+2012
           - Ccg(:,yp,ol,sp).*sth(:,yp,ol,sp).*E(:,yp,ol,sp).*delx(:,yp,ol,sp))...
           ./((dely(:,yp,ol,sp) + dely(:,:,ol,sp)).*(delx(:,yp,ol,sp) + delx(:,:,ol,sp)))*4;
      
       advect(:,:,ol,sm) = advect(:,:,ol,sm)+(Ccg(:,ym,ol,sm).*sth(:,ym,ol,sm).*E(:,ym,ol,sm).*delx(:,ym,ol,sm)...            % advection term in eqn. B6, Donelan+2012
           - Ccg(:,:,ol,sm).*sth(:,:,ol,sm).*E(:,:,ol,sm).*delx(:,:,ol,sm))...
           ./((dely(:,:,ol,sm) + dely(:,ym,ol,sm)).*(delx(:,:,ol,sm) + delx(:,ym,ol,sm)))*4;
      
% -- Full energy from source, sink, and advection -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       E1(:,:,ol,:) = E1(:,:,ol,:) - newdelt*advect(:,:,ol,:);                                                                % eqn. B6, Donelan+2012
      
       % clean up
       E1(D <= 0) = 0;                                                                                                        % energy on land is set to zero
       E1=real(E1);                                                                                                      
       E1(E1 < 0)=0;
      
% -- Compute refraction including wave-current interaction ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      Crot = zeros(size(E));
       Crot(:,:,ol,:) = (c(xm,:,ol,:)-c(xp,:,ol,:)).*sth(:,:,ol,:)./(delx(xm,:,ol,:)+delx(xp,:,ol,:))...                      % eqn. A8, Donelan+2012
           -(c(:,ym,ol,:)-c(:,yp,ol,:)).*cth(:,:,ol,:)./(dely(:,ym,ol,:)+dely(:,yp,ol,:));
       % Determine if rotation is clockwise or counter clockwise
       Crotcw = zeros(size(Crot));Crotccw = zeros(size(Crot));                  
       Crotccw(Crot>0) = Crot(Crot>0); Crotcw(Crot<0) = Crot(Crot<0);
       Crotccw(Crotccw > dth/newdelt) = dth/newdelt;
       Crotcw(Crotcw < -1*dth/newdelt) = -1*dth/newdelt;

       E1(:,:,ol,cw) = E1(:,:,ol,cw) - newdelt*Crotcw(:,:,ol,:).*E(:,:,ol,:)/dth;                                                          % eqn. A1, Donelan+2012 (clockwise rotation)
       E1(:,:,ol,:) = E1(:,:,ol,:) + newdelt*Crotcw(:,:,ol,:).*E(:,:,ol,:)/dth;                                                            % add clockwise rotation
       E1(:,:,ol,ccw) = E1(:,:,ol,ccw) + newdelt*Crotccw(:,:,ol,:).*E(:,:,ol,:)/dth;                                                       % eqn. A1, Donelan+2012 (counterclockwise rotation)
       E1(:,:,ol,:) = E1(:,:,ol,:) - newdelt*Crotccw(:,:,ol,:).*E(:,:,ol,:)/dth;                                                           % add counter-clockwise rotation
       E1(E1 < 0) = 0; E1(D <= 0) = 0;
       E = real(E1);E(isnan(E)) = 0;

       % Make all land points = 0 and upstream values equal to coarse grid
       E(D <= 0) = 0;
                
       % Compute form stress spectrum.
       tauE = sum(wn.*E.*Sin.*cth./c,4)*dth;                                                                                            % eastward wind stress (eqn. 5, Donelan 2012)
       tauN = sum(wn.*E.*Sin.*sth./c,4)*dth;                                                                                            % northward wind stress (eqn. 5, Donelan 2012)
       % Add wind speed dependent tail of slope, "mtail" pinned to the highest wavenumber, "wnh".
       mtail = 0.000112*U_10.*U_10 - 0.01451.*U_10 - 1.0186;
       wnh = squeeze(wn(:,:,model.Fdim,1));                                                                                                % largest wavenumber
      
       tauE = planet.rho_liquid*(sum((planet.gravity+planet.surface_tension.*squeeze(wn(:,:,:,1)).^2./planet.rho_liquid).*squeeze(dwn(:,:,:,1)).*tauE,3) + ...
           (planet.gravity+planet.surface_tension.*squeeze(wn(:,:,model.Fdim,1)).^2./planet.rho_liquid).*squeeze(tauE(:,:,model.Fdim)).*wnh.^(-mtail).*(kutoff.^(mtail+1)-wnh.^(mtail+1))./(mtail+1));    % eqn. 5, Donelan+2012
       tauN = planet.rho_liquid*(sum((planet.gravity+planet.surface_tension.*squeeze(wn(:,:,:,1)).^2./planet.rho_liquid).*squeeze(dwn(:,:,:,1)).*tauN,3) + ...
           (planet.gravity+planet.surface_tension.*squeeze(wn(:,:,model.Fdim,1)).^2./planet.rho_liquid).*squeeze(tauN(:,:,model.Fdim)).*wnh.^(-mtail).*(kutoff.^(mtail+1)-wnh.^(mtail+1))./(mtail+1));    % eqn. 5, Donelan+2012
      
       Cd = abs(tauE + i*tauN)./rhoa./(U_z.^2);                                                                                             % shear stress and law of the wall
       Cdf = Cd;
       Ustar_smooth = smooth_nu((1-wfac)*U_z(:),model.z_data,planet.nua);
       Ustar_smooth = reshape(Ustar_smooth,model.LonDim,model.LatDim);                                                                                % Surface current (friction velocity for hydraulically smooth interface)
      
       Ustar_smooth = (Ustar_smooth./U_z).^2;
       Cds =  Ustar_smooth;
      
       Ustar_smooth = U_z.^2.*(0.3333*Ustar_smooth + 0.6667*(Ustar_smooth.^2)./(Ustar_smooth + Cd));
       tauE = tauE + rhoa*Ustar_smooth.*cos(squeeze(windir(:,:,1,1)));                                                                      % wind stress + wind momentum in Eastward direction
       tauN = tauN + rhoa*Ustar_smooth.*sin(squeeze(windir(:,:,1,1)));                                                                      % wind stress + wind momentum in Northward direction

       Cd = abs(tauE + i*tauN)./rhoa./(U_z.^2);                                                                                             % form drag coefficient (eqn. 8, Donelan+2012)
       Cd = clip(Cd,0,2.5e-3);                                                                                                              % put a cap here on the absolute size of the drag coefficient so U at lambda/2 doesn't go negative and become unphysical, max from Donelan+2004, fig. 2

% -- Sig wave height -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      % integrate spectrum to find significant height  
            ht = sum(dwn.*wn.*E,4)*dth;                            % integral = sum(wavespectrum*wn*del_wn) -->  spectral moment
            ht = sum(ht,3);                                        % sum the prev sum over all frequencies to get the zeroth order moment (aka variance of sea surface (1/2a^2))
            ht = 4*sqrt(abs(ht));                                  % significant wave height (from zeroth order moment of surface)
            ks = ht;
            sigH(t) = ht(model.long,model.lat);                    % return significant wave height at specified lat,lon coordinates 

            if nargout > 1
                htgrid{t} = ht;                                    % return significant wave height at each spatial point (m,n) on the grid
            end
          

% -- mean slope ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
       % integrate spectrum to find mean slope                                                                                                                         
        ms = sum(dwn.*wn.^3.*E,4)*dth;                                                                                                              % slope = angle water surface makes with flat surface
        ms = sum(ms,3);
        ms = sqrt(ms);                                                                                                                              % standard deviation of water surface = sqrt(variance of water surface)

       % Integrate spectrum over wavenumber to plot directional plot of 10th wavelength.
       Wavel = sum(squeeze(wn(:,:,10,:)).*squeeze(E(:,:,10,:).*dwn(:,:,10,:)),3)./sum(squeeze(E(:,:,10,:).*dwn(:,:,10,:)),3);                           % integrating over frequency spectrum for plotting
       Wavel = 2*pi./Wavel;
       Wavel = Wavel.^(0.25); 
       % Direction of one wavelength.
       KD = squeeze(sum(dwn(:,:,10,:).*E(:,:,10,:),3));                                                                                                 % direction vector of wave at 10th frequency bin (near center for o = 25)
       Dir_sin = sum(KD.*squeeze(sin(waveang(:,:,10,:))),3)./sum(KD,3);                                                                                 % x-component of wave at 10th frequency bin
       Dir_cos = sum(KD.*squeeze(cos(waveang(:,:,10,:))),3)./sum(KD,3);                                                                                 % y-component of wave at 10th frequency bin
       mdir = atan2(Dir_sin,Dir_cos);                                                                                                                   % average wave direction [-pi to pi]
      
       cgmax = max(max(max(max((E>1e-320).*Cg))));                                                                                                      % Adjust cgmax for active energy components.
      

   end % end while loop for sub-time steps
  
   fraction_time_completed = t/model.num_time_steps;
   fprintf('u = %.2f: fraction time completed: %.2f\n',UU,fraction_time_completed);

  if nargout > 2
    E_each{t} = E;    % return energy spectrum for each wind speed  at each time step t
  end

  if Etc.savedata
    save([TitanResults,'/New_',int2str(file),'_',int2str(t)],'E','ht','freqs','oa','Cd','Cdf','Cds','Sds','Sds_wc','Sin','Snl','Sdt','Sbf','ms')
  end
   
   
   if t > 100 && sigH(t-1)/sigH(t) < model.tolH                                                                                                          % will break out of wind speed loop if waves haven't changed by more than the tolerance level tolH 
       disp('Waves have reached 99% of maturity.')
       break
   elseif t > 1
       fprintf('t_n-1/t_n: %.6f\n',sigH(t-1)/sigH(t));
       if sigH(t-1)/sigH(t) > 1
            warning('Numerical ringing. Suggested fix: (1) Re-run with a larger cutoff frequency (2) Re-run with a larger maximum time step.')
       end
   end

   

end % end of loop for Tsteps

[~,minind] = min(abs(wind.speed-UU));
disp(['Finished Wind Speed ' num2str(UU) ' m/s (' num2str(minind) ' Of ' num2str(numel(wind.speed)) ' )']);
  

%% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if Etc.savedata
    % Make gif of sig height results:
    idx_end = idx;
    filename = 'SigH.gif';                                                 % Specify the output file name
    for idx = 1:idx_end-1
       [A,map] = rgb2ind(im{idx},256);
       if idx == 1
           imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
       else
           imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
       end
    end
end
disp('================================================================')
disp('Model Run Complete')
disp('================================================================')
toc
end
% ==============================================================================================================================================================================================================================================================================================================================================================================================================================
% User-Defined Functions:   

%% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function ustar = smooth_nu(U,z,nu)
% To calculate the friction velocity for smooth flow of any gas.
%
% function ustar = smooth(U,z)
%
% ustar = friction velocity [m/s]
% U     = wind speed [m/s] at height z
% z     = height of wind speed measurement [m]
% nu    = kinematic viscosity of gas [m^2/s]

    ustar = NaN(1,length(U));
    z0 = 0.001;

    for j = 1:length(U)
       m = 0;    
       for k = 1:6
           m = m+1;
           ust = 0.4*U(j)/(log(z/z0));
           z0 = (1/9)*nu/ust;
       end
       ustar(j) = ust;
    end
end

function TitanResults = make_log(planet,model,wind,uniflow,Etc)
% make a log and save file location for result outputs

    assert(rem(model.Dirdim,8)==0,'Model input parameter p must be factorable by 8.')
    % -- prepare log file for commands 
    dfile=strcat(string(datetime('now','TimeZone','local','Format','ddMMyy_HHmmss')),'_wind_speed_',num2str(wind.speed),'_RunLog.txt');
    diary(dfile);
    RAII.diary = onCleanup(@() diary('off'));                                  % auto-closes logging function on error
    % -- create output directory for results 
    TitanResults = strcat('wind_speed_',num2str(wind.speed)); 
    
    if exist(TitanResults, 'dir') == 7  && Etc.savedata                        % make output directory 'Results' if doesn't already exist
       oldmatfiles = fullfile(TitanResults, '*.mat');                          % empties output directory from previous runs
       oldmatloc = dir(oldmatfiles);
       for kk = 1:length(oldmatloc)
           basemat = oldmatloc(kk).name;
           fullmat = fullfile(TitanResults,basemat);
           fprintf(1,'Deleting previous .mat files %s\n',fullmat);
           delete(fullmat);
       end
    elseif Etc.savedata
	    mkdir(TitanResults);
    end

% adds run details to console and saves to log file 
    disp('================================================================')
    disp(['Directional Wave Spectrum -- last updated: ' dir('makeWaves.m').date])
    disp(['Wind Speed(s) to Run: ' regexprep(mat2str(wind.speed),{'\[', '\]', '\s+'}, {'', '', ','}) ' m/s']);
    disp('================================================================')

    fprintf('Planet:\t\t\t\t %s\n',planet.name)
    fprintf('\trho_liquid:\t\t %.3f\n',planet.rho_liquid)
    fprintf('\tnu_liquid:\t\t %.3e\n',planet.nu_liquid)
    fprintf('\tnua:\t\t\t %.3e\n',planet.nua)
    fprintf('\tgravity:\t\t %.3f\n',planet.gravity)
    fprintf('\tsurface_temp:\t\t %.3f\n',planet.surface_temp)
    fprintf('\tsurface_press:\t\t %.3f\n',planet.surface_press)
    fprintf('\tsurface_tension:\t %.3f\n',planet.surface_tension)
    
    fprintf('\n')
    
    fprintf('Model\n')
    fprintf('\tm:\t\t\t %i\n',model.LonDim)
    fprintf('\tn:\t\t\t %i\n',model.LatDim)
    fprintf('\to:\t\t\t %i\n',model.Fdim)
    fprintf('\tp:\t\t\t %i\n',model.Dirdim)
    fprintf('\tlong:\t\t\t %i\n',model.long)
    fprintf('\tlat:\t\t\t %i\n',model.lat)
    fprintf('\tgridX:\t\t\t %i\n',model.gridX)
    fprintf('\tgridY:\t\t\t %i\n',model.gridY)
    fprintf('\tmindelt:\t\t %.2e\n',model.mindelt)
    fprintf('\tmaxdelt:\t\t %.2f\n',model.maxdelt)
    fprintf('\ttime_step:\t\t %i\n',model.time_step)
    fprintf('\tnum_time_steps:\t\t %i\n',model.num_time_steps)
    fprintf('\ttolH:\t\t\t %.2f\n',model.tolH)
    fprintf('\tcutoff_freq:\t\t %i\n',model.cutoff_freq)
    fprintf('\tmin_freq:\t\t %.2f\n',model.min_freq)
    fprintf('\tmax_freq:\t\t %.2f\n',model.max_freq)
    fprintf('\tz_data:\t\t\t %.2f\n',model.z_data)
    fprintf('\ttune_A1:\t\t %.2f\n',model.tune_A1)
    fprintf('\ttune_mss_fac:\t\t %.2f\n',model.tune_mss_fac)
    fprintf('\ttune_Sdt_fac:\t\t %.3f\n',model.tune_Sdt_fac)
    fprintf('\ttune_Sbf_fac:\t\t %.3f\n',model.tune_Sbf_fac)
    fprintf('\ttune_cotharg:\t\t %.2f\n',model.tune_cotharg)
    fprintf('\ttune_n:\t\t\t %.2f\n',model.tune_n)
    fprintf('\tbathy_map (deepest):\t %.2f\n',max(max(model.bathy_map)))
    
    allSame = all(all(model.bathy_map == model.bathy_map(1,1)) == 1);
    if allSame
        slopex = 0;
        slopey = 0;
    else
        slopex = max(max(model.bathy_map))/(model.LonDim*model.gridX);
        slopey = max(max(model.bathy_map))/(model.LatDim*model.gridY);
    end
    
    
    fprintf('\tbathy_map (slopex):\t %.2e\n',slopex)
    fprintf('\tbathy_map (slopey):\t %.2e\n',slopey)
    
    
    fprintf('\n')
    
    fprintf('Wind\n')
    [fprintf('\tspeed: \t\t\t'),fprintf(' %.2f,', wind.speed(1:end-1)), fprintf(' %.2f\n', wind.speed(end))];
    fprintf('\tdir: \t\t\t %.2f',wind.dir)
    
    fprintf('\n')
    
    fprintf('Uniflow\n')
    fprintf('\tEast:\t\t\t %.2f\n',uniflow.East)
    fprintf('\tNorth:\t\t\t %.2f\n',uniflow.North)
    
    fprintf('\n')
    
    fprintf('Etc\n')
    fprintf('\tshowplots:\t\t %i\n',Etc.showplots)
    fprintf('\tsavedata:\t\t %i\n',Etc.savedata)
    fprintf('\tshowlog:\t\t %i\n',Etc.showlog)

disp('================================================================')

end

function surf_extrema(my_array,variable_name,model,extrema_type)
% make a surf plot of the max and min values of a variable of interest

    %my_array = squeeze(my_array(round(model.LonDim/2),round(model.LatDim/2),:,:));

    reshaped_array = reshape(my_array,model.LonDim*model.LatDim,model.Fdim*model.Dirdim);
    
    if strcmp(extrema_type,'tallest')
        val_ar = max(reshaped_array,[],2);
    elseif strcmp(extrema_type,'smallest')
        val_ar = min(reshaped_array,[],2);
    else
        error('SURF_EXTREMA: extrema type not specified')
    end

    valgrid = reshape(val_ar,model.LonDim,model.LatDim);
    

    figure;
    surf(valgrid,'FaceColor','interp')
    colorbar
    title([variable_name,' ',extrema_type,' at center of grid'])
    view(2)
    hold on

end

function plot_freq_depend(my_array,variable_name,depth,frequencies,model)
% plot the frequency dependence of a variable in the deepest and shallowest section of the grid

    [deep_spot_depth, deep_spot_li] = max(depth(:));
    [deep_ri,deep_ci] = ind2sub(size(depth),deep_spot_li);

    [shallow_spot_depth, shallow_spot_li] = min(depth(:));
    [shallow_ri,shallow_ci] = ind2sub(size(depth),shallow_spot_li);



    figure;
    plot(frequencies,squeeze(my_array(deep_ri,deep_ci,:,round(model.Dirdim/2))),'LineWidth',3)
    hold on
    plot(frequencies,squeeze(my_array(shallow_ri,shallow_ci,:,round(model.Dirdim/2))),'LineWidth',3)
    xline(frequencies(model.cutoff_freq),'-','cutoff frequency')
    legend({strcat('Deep=',num2str(deep_spot_depth)), strcat('Shallow=',num2str(shallow_spot_depth))},'Location', 'Best')
    xlabel('frequency [Hz]')
    ylabel(variable_name)
    set(gca, 'XScale', 'log')
    hold off;
    

end
% ==============================================================================================================================================================================================================================================================
% REFERENCES CITED IN COMMENTS:
% 1. Donelan+2012 : Donelan et al. 2012 "Modeling waves and wind stress" (JGR)
% 2. Donelan+2001 : Donelan et al. 2001 "A Nonlinear Dissipation Function due to Wave Breaking (Proc. ECMWF Workshop on Ocean Wave Forecasting)
% 3. Kinsman      : Kinsman, Blair 1965 "Wind Waves: Their Generation and Propagation on the Ocean Surface"
% 4. Donelan+2004 : Donelan et al. 2004 "Aerodynamic Roughness of the Ocean in Strong Winds"
% ==============================================================================================================================================================================================================================================================