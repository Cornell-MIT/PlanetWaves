function [sigH,T0] = makeWaves(windspeeds,rho_liquid,nu_liquid,bathy_map,time_step_size,num_time_steps)
%% ==========================================================================================================================================================================================================================================================================
%% ==========================================================================================================================================================================================================================================================================
% This code calculates E(x,y,k,theta) for wave field using an energy balance between wind-input and multiple dissipation terms (see Donelan et al. 2012 Modeling Waves and Wind Stress).
%
% The equation to solve is:
%   E_{n+1} = E_{n} + del*(-Cg*cos(theta)*dE/dx - Cg*sin(theta)*dE/dy] + Sin + Snl - Sds
%

% Dimensions of x,y,k,th are m,n,o,p
%
%   Arguments:
%       windspeeds     : speed of wind [m/s] above boundary layer
%       rho_liquid     : density of the liquid [kg/m3]
%       nu_liquid      : kinematic viscocity of liquid [m2/s]
%       bathy_map      : m x n array of depth [m] (+ values = subsurface, - values = subaerial)
%       time_step_size : maximum size of time step
%       num_time_steps : total number of time steps to iterate over 
%
%   Returns:
%       sigH           : largest signifigant wave height over the entire spatial array [m]
%       dwn            : dominant wavenumber [1/m]
%   
%       
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% External function requirements:
%   (1) none
%
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Aprox time to run: 20 minutes (for 2 wind speeds, Time = 100, Tsteps = 200)
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Cite:
% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% Authors: Mark Donelan, Alex Hayes, Charlie Detelich, Una Schneck
%
%% ==========================================================================================================================================================================================================================================================================
%% ==========================================================================================================================================================================================================================================================================
% -- prepare log file for commands -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dfile=[datestr(now,'mmddyyyy_HHMMSS'),'_FetchLaws.txt'];
diary(dfile);
diary on;
% -- create output directory for results -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
parentFolder = fileparts(pwd);
TitanResults = sprintf('%s\\Titan%s', pwd);
if ~exist(TitanResults, 'dir')
 mkdir(TitanResults);
end
% -- prepare .mat files to be saved during loops -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
file = 1;
filec = 1;
% -- SHOW PLOTS  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
showplots = 1;                                                             % 0 = no plots made, 1 = plots every tenth loops
% -- MODEL SET UP ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% time inputs for model
numdays = 1;
Time = time_step_size;                                                     % Total size of time step (maximum size of dynamic sub time step)
Tsteps = num_time_steps;                                                   % maximum number of time steps
mindelt = 0.01;                                                            % Minimum time step to minimize oscillations of delt to very small values.
maxdelt = 5000.0;%2000.0;                                                          % Maximum time step to prevent large values as wind ramps up at startup
tolH = 0.001;                                                              % minimum change in SigH to stop model at (otherwise runs to full Tsteps)
% global constants
kappa = 0.4;                                                               % Von-Karman constant
i = sqrt(-1);                                                              % imaginary number i
kgmolwt = 0.028;                                                           % gram molecular weight [Kgm/mol]
RRR = 8.314;                                                               % Universal gas constant [J/K/mol]
% spatial grid
[m,n] = size(bathy_map);                                                   % [number of gridpoints in y-direction,number of gridpoints in y-direction]
% frequency and direction bins
o = 25;                                                                    % number of frequency bins
p = 64;                                                                    % number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
dr = pi/180;                                                               % conversion from degrees to radians
dthd = 360/(p);                                                            % small angle for integration (for angular direction) [degrees]
dth = 360/p;                                                               % small angle for integration [degrees]
% Set up geographic deltas
Delx = 1000.0;                                                             % grid step size [m]
Dely = 1000.0;                                                             % grid step size [m]
dely = Dely*ones(m,n,o,p);
delx = Delx*ones(m,n,o,p);
% Wind:
gust = 0;                                                                  % Gust factor applied to wind at each time step
zref = 20;                                                                 % height of reference wind speed, normally at 20m [m]
z = 10;                                                                    % height of measured wind speed [m]
wfac = 0.035;                                                              % winddrift fraction of Uz ffor U10m
UUvec = windspeeds;                                                        % wind speeds to model [m/s]
wind_angle = 0;                                                            % wind approach angle
% currents:
uec = 0;                                                                   % Eastward current [m/s]
unc = 0;                                                                   % Northward current [m/s]
% location of interest within grid
long = 10;                                                                 % longitude within grid for Punga Mare
lati = 10;                                                                 % latitude within grid for Punga Mare
% Surface Conditions:
sfcpress = 1*101300;%1.5*101300;                                                     % surface pressure [Pa]
sfctemp = 273;%92;                                                              % surface temperature [K]
sfcT = 0.072;%0.018;                                                              % surface tension of liquid [N/m] (0.027 for 75% methane)
% Ideal gas law: PV = nRT
% Densities:
rhoa = sfcpress*kgmolwt/(RRR*sfctemp);                                     % air density [kg/m3]
rhow = rho_liquid;                                                         % liquid density [kg/m3] (550 kg/m^3 for 75% methane)
rhorat=rhoa/rhow;                                                          % air-water density ratio.
g = 9.81;%1.35;                                                                  % gravity [m/s2]
% Viscocities:
nu = nu_liquid;                                                            % liquid viscocity [m2/s]
nua = 0.0126/1e4;                                                          % atmospheric gas viscosity [m2/s]
% Wavenumber limits:
kutoff = 1000;                                                             % wavenumber cut-off due to Kelvin-Hemholtz instabilities
kcg = sqrt(g*rhow/sfcT);                                                   % wavenumber of slowest waves defined by the capillary-gravity waves, from Airy dispersion: omega^2 =gktanh{k)
% modified wavenumbers
kcgn = 1.0*kcg;                                                            % n power min is not shifted with respect to kcg
kcga = 1.15*kcg;                                                           % a min is shifted above kcg.
% frequency limits:
f1 = 0.05;                                                                 % minimum frequency
f2 = 35;                                                                   % maximum frequency
% create frequency limits for spectrum
dlnf=(log(f2)-log(f1))/(o-1);                                              % frequency range
f = exp(log(f1)+[0:o-1]*dlnf);                                             % frequencies for spectrum
dom = 2*pi*dlnf.*f;                                                        % dominant angular frequency (w = 2pi*f)
freqs = f;                                                                 % save a copy of frequencies
% Frequency bins:
o = length(f);                                                             % set o to maximum frequency
oa = 7;                                                                    % top bin advected (make wavelengths [vector wn] at oa 1/20th of the grid spacing)
ol = [1:oa];                                                               % bins for long frequencies
os = [oa+1:o];                                                             % bins for short frequencies
% Compute diffusion values in 2 freqs and 2 directions:
%   bf1 + bf2 = 1, and bt1 + bt2 = 1.
bf1 = exp(-15.73*dlnf*dlnf);                                               % part of non-linear source term for downshifting and spilling breakers, eqn. 21 Donelan 2012
bf2 = exp(-15.73*4*dlnf*dlnf);                                             % part of non-linear source term for downshifting and spilling breakers, eqn. 21 Donelan 2012
bf1a = bf1/(bf1 + bf2);                                                    % normalization (eqn. 21, Donelan 2012)
bf2 = bf2/(bf1 + bf2);                                                     % normalization (eqn. 21, Donelan 2012
bf1 = bf1a;
% Compute Snl_fac
A = exp(dlnf);
B = A.*A;
fac = bf1.*(1-1/B) + bf2.*(1-1./B.^2);                                     % eqn. 21, Donelan 2012
Snl_fac = ((1.0942/A)^1.9757)/fac;                                         % eqn. 21, Donelan 2012

waveang = ([0:p-1]-p/2+0.5)*dthd*dr;                                       % wave angle [radians]
op = [p/2+1:p 1:p/2];                                                      % opposite directions from the p array for the computation of MSS colinear (with and against) waves.
orm = [p:-1:1];                                                            % reflected directions for 100% reflection from y-boundaries.
th = ([0:p-1]-p/2+0.5)*dthd*dr;                                            % angle of wave propogation (phi in Donelan 2012) 
cth=cos(th);                                                               % cosine of angle in wave propogation direction 
sth=sin(th);                                                               % sine of angle in wave propogation direction
dth=dthd*dr;                                                               % 2pi/p (small angle for integration [radians])
% compute cos^2 for calculation of mss vs direction.
tth = ([0:p-1])*dthd*dr;                                                   % angular difference between short waves and longer waves
cth2 = cos(tth).^2;                                                        % cosine square of angular difference between short waves and longer waves
% indices for refraction rotation:
cw = ([p 1:p-1]);                                                          % clockwise indices
ccw = ([2:p 1]);                                                           % counterclockwise indices
% upwave indices for advective term:
xp = [1 1:m-1];
yp = [1 1:n-1];
xm = [2:m m];
ym = [2:n n];
% Index being used for advection:
cp = find(cth > 0);
cm = find(cth < 0);
sp = find(sth > 0);
sm = find(sth < 0);
% reshape matrices
cth = repmat(cth,[o 1 m n]);
cth=shiftdim(cth,2);
sth = repmat(sth,[o 1 m n]);
sth=shiftdim(sth,2);
cth2 = repmat(cth2,[o 1 m n]);
cth2=shiftdim(cth2,2);
%% -- Lake Geometry ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

D = bathy_map;                                                             % depth of liquid [m]

j12 = find(D<=0);                                                          % turn all negative depths to zero
D(j12) = 0;
% Impose depth limits
jj = find (D < 0);
D(jj) = 0;                                                                 % limit land elevations to 0 to avoid dD/dx, dD/dy errors in refraction calculation
D(:,1) = 0; D(1,:) = 0; D(:,end) = 0; D(end,:) = 0;                        % Set array boundary depths to 0 (absorbtive boundaries)
% plot the bathymetry
[xplot,yplot] = meshgrid(1:m,1:n);
if showplots
    figure;
    surf(xplot,yplot,D','EdgeColor','k')
    myc = colorbar;
    myc.Label.String = 'Liquid Depth [m]';
    title('Lake Model Bathymetry')
end
%% -- Fractions of terms -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Fractions applied to terms from fitting with experiment/observations:
mss_fac = 240;                                                             % MSS adjustment to Sdiss. 2000 produces very satisfactory overshoot, 400 does not. (A4 in Donelan 2012; eqn. 6 Donelan 2001)
% Snl_fac = 3.5;                                                           % fraction of Sdiss that goes to longer wavenumbers
% Snl_fac = 4.3;                                                           % fraction of Sdiss that goes to longer wavenumbers
Sds_fac = 1.0;                                                             % fraction with variable power nnn in Sds that goes into spectrum
Sdt_fac = 0.001;                                                           % fraction of Sdt that goes into spectrum (A4 in eqn. 20, Donelan 2012; A4 = 0.01?)
Sbf_fac = 0.002;                                                           % fraction of Sbf that goes into spectrum (0.004 is for smooth sandy bottom) (Gf in eqn. 22, Donelan 2012) 
%% -- wavemumber and Power of Sds ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
wn(:,:,:) = ones(m,n,o);
nnn(:,:,:) = ones(m,n,o);
ann(:,:,:) = ones(m,n,o);

for jm = 1:m
   for jn = 1:n
       if D(jm,jn) > 0
           wn(jm,jn,:) = wavekgt(f,D(jm,jn),g,sfcT,rhow,1e-4);                                                                 % wave number (using linear wave dispersion)
           nnn(jm,jn,:) = 1.2 + 1.3*(abs(2 - (1+3*(wn(jm,jn,:)./kcgn).^2)./(1+(wn(jm,jn,:)./kcgn).^2)).^2.0);                  % Power n of Sds on the degree of saturation [eqn. 15, Donelan 2012, eqn. 5, Donelan 2001]
           ann(jm,jn,:) = 0.04 + 41.96*(abs(2 - (1+3*(wn(jm,jn,:)./kcga).^2)./(1+(wn(jm,jn,:)./kcga).^2)).^4.0);               % Power of Sds
       end
   end
end
% reshape the matrix
wn = repmat(wn,[1 1 1 p]);
nnn = repmat(nnn,[1 1 1 p]);
ann = repmat(ann,[1 1 1 p]);
% conversions for nnn?
nnn = (2.53/2.5).*nnn;
nnninv = 1./nnn;
% reshape the matrix 
Depth = D; D = repmat(D,[1 1 o p]);
f=repmat(f',[1 p m n]);
f=shiftdim(f,2);
dom=repmat(dom',[1 p m n]);
dom=shiftdim(dom,2);
%% -- ocean currents ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Uer=uec*D + 0.0;                                                                                                                           % Eastward current, m/s
Uei=unc*D + 0.0;                                                                                                                           % Northward current, m/s
%% -- wave speed and group velocity --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
c=(2*pi*f)./wn;                                                                                                                           % phase speed
c(D<=0) = 0;                                                                                                                              % phase speed on land is set to zero
Cg=zeros(size(c));                                                                                                                        % initialize group speed
Cg(D>0) = c(D>0)./2.*(1 + 2*wn(D>0).*D(D>0)./sinh(2*wn(D>0).*D(D>0)) + 2*sfcT.*wn(D>0)./rhow./(g./wn(D>0) + sfcT.*wn(D>0)./rhow));        % Group velocity for all waves (Kinsman "Wind Waves: Their Generation and Propogation on the Ocean Surface")
dwn = ones(m,n,o,p);                                                                                                                      % initalize dominant wavenumber (c = dw/dk)
dwn(D>0)=dom(D>0)./abs(Cg(D>0));                                                                                                          % remove any values on land for dominant wavenumber (c = dw/dk)
Cg(D<=0) = 0;                                                                                                                             % zero all the group velocities on land
l2=abs(c)./f/2;                                                                                                                           % wavelength/2
lz=abs(c)./f/2/pi;                                                                                                                        % wavelength/2/pi: kz = 1 for drift current action
% Set l2, lz > zref equal to zref
l2(l2 > zref) = zref;
lz(lz > zref) = zref;
E = zeros(m,n,o,p);                                                                                                                       % initializing the energy term in x,y,k,theta with zeros
Ccg = Cg + Uer.*cth + Uei.*sth;                                                                                                           % total wave group velocity = wave group velocity + currents
cgmax = max(max(max(max(Cg))));                                                                                                           % fastest group velocity
mindelx = min(squeeze(delx(1,:,1,1)));                                                                                                    % smallest spatial grid spacing
delt = 0.7 * mindelx / cgmax;                                                                                                             % advection limited Courant condition.
waveang = repmat(waveang,[o 1 m n]);
waveang=shiftdim(waveang,2);
eval(['save Titan/Titan_',int2str(file),'_ref m n o p freqs f wn dwn D Uei Uer Cg c delx dely dthd waveang'])
file = file - 1;
%% -- loop through wind speeds --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tic
sumt = 0;                                                                  % initialize total model time passed
tplot = -1;
idx = 1;                                                                   % for frame for gif
sigH = NaN(numel(UUvec),length(1:Tsteps));                                 % initialize sigH for returning 

for iii=1:numel(UUvec)                                                     % loop through wind velocities
   UU = UUvec(iii);                                                        % UU = wind speed for loop
   file = file + 1;                                                        % for naming files
  
   U = UU*ones(m,n);                                                       % set wind velocity everywhere in x-y plane
   windir = wind_angle*ones(m,n);                                          % set wind direction everywhere in x-y plane
   windir = repmat(windir,[1 1 o p]);                                      % reshape wind direction matrix
  
  
   U_z = U + 0.005;                                                        % wind at modelled height (plus small number to wind speed to avoid division by zero)
  
   % drag coefficient 
   Cd = 1.2*ones(size(U));                                                 % drag coefficient for weak/moderate winds
   uj = find(U > 11);
   Cd(uj) = 0.49 +0.065*U(uj);                                             % drag coefficient for strong winds
   Cd=Cd/1000;
  
   modt = 0;                                                               % model time initializiation
  
  
   % Currents superimposed onto wave-driven flow
   Uer=0*D + 0.0;                                                          % Eastward current [m/s]
   Uei=0*D + 0.0;                                                          % Northward current [m/s]
  
%% -- loop through model time to grow waves --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   for t = 1:Tsteps                                                                                                                                                 % loop through time
      
       Newdelt = [];  
       sumt = 0;                                                                                                                                                    % intitalize total time within timestep t
       tplot = - 1; 
       while ((Time - sumt) > 0)                                                                                                                                    % each time step is determined as min[max(0.1/(Sin-Sds) 0.0001) 2000 UserDefinedTime CourantGrid], iterate until the user-defined time is larger than the time passed within timestep t
          
           if gust > 0
                U = U.*(1+gust*randn(size(U)));                                                                                                                     % add random bursts of gusts to wind speed
           end
           tplot = tplot + 1;
           explim = 0.1;                                                                                                                                            % limit for making dynamic time step if the source = dissipation
           delt = 0.7 * mindelx / cgmax;                                                                                                                            % advection limited Courant condition.
          
           Ul = real(U_z.*sqrt(Cd));                                        
           Ustar = Ul;
           U_10 = 2.5*Ul.*log(10/z) + U_z;                                                                                                                          % law of wall
           ustarw = Ustar .* sqrt(rhorat);                                                                                                                          % shear stress
          
           Ul = repmat(Ul,[1 1 o p]);
           U = repmat(U_z,[1 1 o p]);
           Ul = 2.5*Ul.*log(l2/z)+U;
           clear Cd
          
           ustw = repmat(ustarw,[1 1 o p]);
           Ud = - ustw./kappa.*log(lz/z);                                                                                                                           % depth averaged air velocity (law of the wall)
          
          
           % calculate wind input as a fraction of stress
           Ul(D>0) = abs(Ul(D>0).*cos(windir(D>0)-waveang(D>0))-c(D>0)-Ud(D>0).*cos(windir(D>0)-waveang(D>0))-Uer(D>0).*cth(D>0)-Uei(D>0).*sth(D>0))...             % Sin = A1*Ul*((k*wn)/g)*(rhoa/rhow)*E  [eqn. 4, Donelan 2012]
               .*(Ul(D>0).*cos(windir(D>0)-waveang(D>0))-c(D>0)-Ud(D>0).*cos(windir(D>0)-waveang(D>0))-Uer(D>0).*cth(D>0)-Uei(D>0).*sth(D>0));
          
          
           % Adjust input to lower value when waves overrun the wind because as wave speed approaches wind speed, energy extraction becomes less efficient
           Ula = Ul;
           Ul = 0.11*Ul;
           fij=find(Ul<0);
           Ul(fij) = Ula(fij)*0.03;                                                                                                                                 % Field scale (Donelan et al, 2012). 0.135 in Lab
           clear fij;                                                     
           fij = find(cos(windir-waveang)<0);                            
           Ul(fij) = Ula(fij)*0.015;                                                                                                                                % Waves against wind
           clear fij;
% -- Sin --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           % input to energy spectrum from wind
           Sin = zeros(size(Ul));
           Sin(D>0)=rhorat*(Ul(D>0).*2*pi.*f(D>0).*wn(D>0)./(g+sfcT.*wn(D>0).*wn(D>0)./rhow));          % eqn. 4, Donelan 2012
           % limits energy going into spectrum once waves outrun wind
           Heavyside = ones(size(E));
           Heavyside(Sin < 0 & E < 1e-320) = 0;
           Sin = Heavyside.*Sin;                                                                       % Heavyside function kills off Sin for negative Sin (energy and momentum transferring from waves to wind) and for negative/zero energy
           clear Ul Heavyside Ula
          
          
           % Set Sin = 0 on land boundaries to avoid newdelt --> 0.
           Sin(D<=0) = 0;
         
           for tj = 1:p
               short(:,:,:,tj) = sum(E.*cth2(:,:,:,rem([1:p]-tj+p,p)+1),4)*dth;                        % energy in each angular bin get the mean square slope (eqn. 16, Donelan 2012)
           end
           short = (cumsum((wn.^3.*short.*dwn),3)-wn.^3.*short.*dwn);                                  % sqrt of mean square slope (eqn. 16, Donelan 2012)
          
% -- Sdt ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           % Calculate turbulent dissipation: Sdt = 4 nu_t k^2.
           Sdt(:,:,:,:) = repmat(Ustar,[1 1 o p]);
           clear Ustar
           Sdt = Sdt_fac*sqrt(rhorat).*Sdt.*wn;                                                                                % eqn 20, Donelan 2012
%-- Sbf -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           Sbf = zeros(m,n,o,p);
           Sbf(D>0) = Sbf_fac*wn(D>0)./sinh(2*wn(D>0).*D(D>0));                                                                % eqn. 22, Donelan 2012
% -- Sds ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           Sds = zeros(size(Sin));
           Sds(D>0)=abs(ann(D>0)*2*pi.*f(D>0).*(1+mss_fac*short(D>0)).^2.*(wn(D>0).^4.*E(D>0)).^(nnn(D>0)));                                    % LH p 712. vol 1 [eqn. 17, Donelan 2012]
          
          
           % Set Sds = 0 on land boundaries to avoid newdelt --> 0.
           Sds(D<=0) = 0;
          
           % Spread Snl to 2 next longer wavenumbers exponentially decaying as distance from donating wavenumber.
           Snl = zeros(m,n,o,p);
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
           Sds(D>0) = coth(0.2*wn(D>0).*D(D>0)).*Sds(D>0) + Sdt(D>0) + Sbf(D>0) + 4*nu*wn(D>0).^2;                               % Add viscous, turbulent and plunging dissipation after calculation of Snl
           % aa = input - dissipation
           aa = Sin(:,:,ol,:) - Sds(:,:,ol,:);                                                                                   % ol = long wavelength waves
           aa = aa(D(:,:,ol,:)>0);                                                                                                    
           aa = max(abs(aa(:)));
           aaa = explim;
          
           if isnan(aa)                                                                                                          % if source = dissipation then denominator for new possible time step is 5e-5
               aa = aaa/maxdelt;
           end
           newdelt = max([aaa/aa mindelt]);                                                                                      % newdelt = delt to give max of 50% growth or decay
           newdelt = min([newdelt maxdelt (Time-sumt) delt]);                                                                    % min[max(0.1/(Sin-Sds) 0.0001) 2000 TotalModelTime CourantGrid]
           fprintf('newdelt: %.2f\n',newdelt);
          
           % add to model time
           sumt = sumt + newdelt;
           modt = modt + newdelt;
          
           fac_exp=(Sin(:,:,ol,:) - Sds(:,:,ol,:));
          
% -- Wave Energy ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           E1 = zeros(size(E));                                                                                               % E^(n+1) in forward Euler differencing)
           E1(:,:,ol,:) = E(:,:,ol,:).*exp(newdelt*fac_exp);                                                                  % Long waves.
          
          
           E1(:,:,ol,:) = E1(:,:,ol,:) + newdelt*Snl(:,:,ol,:);                                                               % eqn. B3, Donelan 2012 (time discretizaton solition for variance spectrum at next time step)

           cath = ones(size(Sds));                                                                                            % horizontal-to-vertical orbital velocity enhancement which leads to more rapid dissipation in shoaling waves relative to deep water spilling breakers
           cath(D>0)=coth(0.2*wn(D>0).*D(D>0));                                                                               % limits the breaker height to depth of shoaling wave ratio [eqn. 17, Donelan 2012 (but A2 = 42 not 0.2?)  
           
           E2(:,:,:,:) = zeros(m,n,o,p);
           fij = find((Sin(:,:,:,:) - 4*nu*wn(:,:,:,:).^2 - Sdt(:,:,:,:) - Sbf(:,:,:,:)) > 0);                               % finds terms where input > dissipation
           E2(fij) = wn(fij).^(-4).*((Sin(fij)-4*nu*wn(fij).^2 - Sdt(fij) - Sbf(fij))./(cath(fij) ...
               .*ann(fij)*2*pi.*f(fij).*(1+mss_fac*short(fij)).^2)).^(nnninv(fij));                                           % LH p.712, vol 1
           E2(D<=0) = 0; E1(D<=0) = 0; 
           E2(E2<0) = 0; E1(E1<0) = 0;
          
           E1(:,:,os,:) = E2(:,:,os,:);                                                                                      % E1 = E2 at small wavelengths (os = short frequency)                                                                                      
          
           clear E2
          
          
           % max_mss = max(max(max(max(short))))                                                                             % max(mss) should be < 0.01

           E1=real(E1);
           ji=find(E1 < 0);
           E1(ji)=0;
           E1(D <= 0)=0;
           E(D <= 0)=0;
          
           E(:,:,os,:) = E1(:,:,os,:);
           E(:,:,ol,:) = (E(:,:,ol,:)+E1(:,:,ol,:))/2;                                                                      % E = (E+E1)/2 (time-splitting for stable integration) [eqn. B5, Donelan 2012]
          
% -- Compute advection term -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           advect=zeros(m,n,o,p);
           Ccg = Cg + Uer.*cth + Uei.*sth;
          
          
           % Upwave advection with account taken of varying delx and dely
           advect(:,:,ol,cp)=advect(:,:,ol,cp)+(Ccg(:,:,ol,cp).*cth(:,:,ol,cp).*E(:,:,ol,cp).*dely(:,:,ol,cp)...                % advection term in eqn. B6, Donelan 2012
               - Ccg(xp,:,ol,cp).*cth(xp,:,ol,cp).*E(xp,:,ol,cp).*dely(xp,:,ol,cp))...
               ./((dely(xp,:,ol,cp) + dely(:,:,ol,cp)).*(delx(xp,:,ol,cp) + delx(:,:,ol,cp)))*4;
          
           advect(:,:,ol,cm)=advect(:,:,ol,cm)+(Ccg(xm,:,ol,cm).*cth(xm,:,ol,cm).*E(xm,:,ol,cm).*dely(xm,:,ol,cm)...            % advection term in eqn. B6, Donelan 2012
               - Ccg(:,:,ol,cm).*cth(:,:,ol,cm).*E(:,:,ol,cm).*dely(:,:,ol,cm))...
               ./((dely(:,:,ol,cm) + dely(xm,:,ol,cm)).*(delx(:,:,ol,cm) + delx(xm,:,ol,cm)))*4;
          
           advect(:,:,ol,sp)=advect(:,:,ol,sp)+(Ccg(:,:,ol,sp).*sth(:,:,ol,sp).*E(:,:,ol,sp).*delx(:,:,ol,sp)...                % advection term in eqn. B6, Donelan 2012
               - Ccg(:,yp,ol,sp).*sth(:,yp,ol,sp).*E(:,yp,ol,sp).*delx(:,yp,ol,sp))...
               ./((dely(:,yp,ol,sp) + dely(:,:,ol,sp)).*(delx(:,yp,ol,sp) + delx(:,:,ol,sp)))*4;
          
           advect(:,:,ol,sm)=advect(:,:,ol,sm)+(Ccg(:,ym,ol,sm).*sth(:,ym,ol,sm).*E(:,ym,ol,sm).*delx(:,ym,ol,sm)...            % advection term in eqn. B6, Donelan 2012
               - Ccg(:,:,ol,sm).*sth(:,:,ol,sm).*E(:,:,ol,sm).*delx(:,:,ol,sm))...
               ./((dely(:,:,ol,sm) + dely(:,ym,ol,sm)).*(delx(:,:,ol,sm) + delx(:,ym,ol,sm)))*4;
          
% -- Full energy from source, sink, and advection -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           E1(:,:,ol,:) = E1(:,:,ol,:) - newdelt*advect(:,:,ol,:);                                                              % eqn. B6, Donelan 2012
          
           % clean up
           E1(D <= 0) = 0;                                                                                                      % energy on land is set to zero
           E1=real(E1);                                                                                                      
           ji=find(E1 < 0);
           E1(ji)=0;
          
% -- Compute refraction including wave-current interaction ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
          
           
           Crot(:,:,ol,:)=(c(xm,:,ol,:)-c(xp,:,ol,:)).*sth(:,:,ol,:)./(delx(xm,:,ol,:)+delx(xp,:,ol,:))...                                      % eqn. A8, Donelan 2012
               -(c(:,ym,ol,:)-c(:,yp,ol,:)).*cth(:,:,ol,:)./(dely(:,ym,ol,:)+dely(:,yp,ol,:));
          
           % Determine if rotation is clockwise or counter clockwise
           Crotcw = zeros(size(Crot));Crotccw = zeros(size(Crot));                  
           Crotccw(Crot>0) = Crot(Crot>0); Crotcw(Crot<0) = Crot(Crot<0);
           jj = find (Crotccw > dth/newdelt);
           Crotccw(jj) = dth/newdelt;
           jj = find (Crotcw < -1*dth/newdelt);
           Crotcw(jj) = -1*dth/newdelt;

           E1(:,:,ol,cw) = E1(:,:,ol,cw) - newdelt*Crotcw(:,:,ol,:).*E(:,:,ol,:)/dth;                                                          % eqn. A1, Donelan 2012 (clockwise rotation)
           E1(:,:,ol,:) = E1(:,:,ol,:) + newdelt*Crotcw(:,:,ol,:).*E(:,:,ol,:)/dth;                                                            % not sure
           E1(:,:,ol,ccw) = E1(:,:,ol,ccw) + newdelt*Crotccw(:,:,ol,:).*E(:,:,ol,:)/dth;                                                       % eqn. A1, Donelan 2012 (counterclockwise rotation)
           E1(:,:,ol,:) = E1(:,:,ol,:) - newdelt*Crotccw(:,:,ol,:).*E(:,:,ol,:)/dth;                                                           % not sure
           clear Crotccw Crotcw Crot
                    
           ji=find(E1 < 0);
           E1(ji)=0;
           E1(D <= 0)=0;
           E=real(E1);
           E(isnan(E))=0;
           E1 = E;
          
          
           % Make all land points = 0 and upstream values equal to coarse grid
           E(D <= 0)=0;
          
           % Make boundaries reflective for tanks.
           % E(:,1,:,:)= E(:,2,:,orm);
           % E(:,n,:,:)= E(:,n-1,:,orm);
           % E(1,:,:,:)= E(2,:,:,orm);
           % E(m,:,:,:)= E(m-1,:,:,orm);
           % E(1,:,:,:)=E(1,:,:,:)+paddle(1,:,:,:);
           % E(D<=0)=0;% Land values.
          
           % Compute form stress spectrum.
           tauE = sum(wn.*E.*Sin.*cth./c,4)*dthd*dr;                                                                                            % eastward wind stress (eqn. 5, Donelan 2012)
           tauN = sum(wn.*E.*Sin.*sth./c,4)*dthd*dr;                                                                                            % northward wind stress (eqn. 5, Donelan 2012)
           % Add wind speed dependent tail of slope, "mtail" pinned to the highest wavenumber, "wnh".
           mtail = 0.000112*U_10.*U_10 - 0.01451.*U_10 - 1.0186;
           wnh = squeeze(wn(:,:,o,1));                                                                                                          % largest wavenumber
          
           tauE = rhow*(sum((g+sfcT.*squeeze(wn(:,:,:,1)).^2./rhow).*squeeze(dwn(:,:,:,1)).*tauE,3) + ...
               (g+sfcT.*squeeze(wn(:,:,o,1)).^2./rhow).*squeeze(tauE(:,:,o)).*wnh.^(-mtail).*(kutoff.^(mtail+1)-wnh.^(mtail+1))./(mtail+1));    % eqn. 5, Donelan 2012
           tauN = rhow*(sum((g+sfcT.*squeeze(wn(:,:,:,1)).^2./rhow).*squeeze(dwn(:,:,:,1)).*tauN,3) + ...
               (g+sfcT.*squeeze(wn(:,:,o,1)).^2./rhow).*squeeze(tauN(:,:,o)).*wnh.^(-mtail).*(kutoff.^(mtail+1)-wnh.^(mtail+1))./(mtail+1));    % eqn. 5, Donelan 2012
           % tauE = rhow*(sum((g+sfcT.*squeeze(wn(:,:,:,1)).^2./rhow).*squeeze(dwn(:,:,:,1)).*tauE,3));
           % tauN = rhow*(sum((g+sfcT.*squeeze(wn(:,:,:,1)).^2./rhow).*squeeze(dwn(:,:,:,1)).*tauN,3));
          
           Cd = abs(tauE + i*tauN)./rhoa./(U_z.^2);                                                                                             % shear stress and law of the wall
           Cdf = Cd;
           Ustar_smooth = smooth_nu((1-wfac)*U_z(:),z,nua);
           Ustar_smooth = reshape(Ustar_smooth,m,n);                                                                                            % Surface current (friction velocity for hydraulically smooth interface)
          
           Ustar_smooth = (Ustar_smooth./U_z).^2;
           Cds =  Ustar_smooth;
          
           Ustar_smooth = U_z.^2.*(0.3333*Ustar_smooth + 0.6667*(Ustar_smooth.^2)./(Ustar_smooth + Cd));
           tauE = tauE + rhoa*Ustar_smooth.*cos(squeeze(windir(:,:,1,1)));                                                                      % wind stress + wind momentum in Eastward direction
           tauN = tauN + rhoa*Ustar_smooth.*sin(squeeze(windir(:,:,1,1)));                                                                      % wind stress + wind momentum in Northward direction
           clear Ustar_smooth;
           Cd = abs(tauE + i*tauN)./rhoa./(U_z.^2);                                                                                             % form drag coefficient (eqn. 8, Donelan 2012)
          
           
           if showplots && rem(tplot,10) == 0                                                                                                   % plot every 10th time step if showplots = 1
              
               close all;

               % diagnostic plot
               figure(16);clf;semilogx(squeeze(wn(2,lati,:,p/2)),squeeze(sum(wn([2:4:m],lati,:,:).*E([2:4:m],lati,:,:),4)*dthd*dr)','*-');grid on
               title(['Omni directional wavenumber spectra along latitude ',num2str(lati),'. Time = ',num2str(modt/3600),' Hours'])
               pause(2)
              
               figure(202);clf;subplot(321);plot(Cd(:,lati).*(U_z(:,lati)./U_10(:,lati)).^2,'.-');% Drag coefficient vs fetch
               title('Drag coefficient')
               grid on
               %pause(2)
               figure(202);hold on;subplot(322);semilogx(squeeze(wn(long,lati,:,10)),squeeze(sum(wn(long,lati,:,:).^2.*Sin(long,lati,:,:).*E(long,lati,:,:),4))*dthd*dr,'*-');
               title('k*Sin');grid on
               %pause(3)
               figure(202);hold on;subplot(323);semilogx(squeeze(wn(long,lati,:,10)),squeeze(sum(wn(long,lati,:,:).^2.*Sds(long,lati,:,:).*E(long,lati,:,:),4))*dthd*dr,'*-');
               figure(202);hold on;subplot(323);semilogx(squeeze(wn(long,lati,:,10)),squeeze(sum(wn(long,lati,:,:).^2.*Sdt(long,lati,:,:).*E(long,lati,:,:),4))*dthd*dr,'*-g');
               figure(202);hold on;subplot(323);semilogx(squeeze(wn(long,lati,:,10)),squeeze(sum(wn(long,lati,:,:).^2.*Sbf(long,lati,:,:).*E(long,lati,:,:),4))*dthd*dr,'*-r');
               title('k*Sds(b), k*Sdt(g), k*Sbf(r)');grid on
               %pause(3)
               figure(202);hold on;subplot(324);semilogx(squeeze(wn(long,lati,:,10)),squeeze(sum(wn(long,lati,:,:).^2.*Snl(long,lati,:,:),4))*dthd*dr,'*-');
               % hold on;subplot(324);semilogx(squeeze(wn(long,lati,:,10)),30*cumsum(squeeze(dwn(long,lati,:,10)).*squeeze(sum(dth*wn(long,lati,:,:).*(Snl(long,lati,:,:)-4*Snl_fac*Sds(long,lati,:,:).*E(long,lati,:,:)),4))),'.r-');
               title('k*Snl');grid on
               % figure(202);hold on;subplot(325);semilogx(squeeze(wn(long,lati,:,10)),squeeze(shel(long,lati,:,1)),'*-');%case 2
               % title('sheltering coeff.');grid on
               figure(202);hold on;subplot(325);semilogx(squeeze(wn(long,lati,:,10)),squeeze(sum(wn(long,lati,:,:).^2.*E(long,lati,:,:),4))*dthd*dr,'*-');
               title('k*spectrum');grid on
           end

% -- Sig wave height -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
          % integrate spectrum to find significant height  
            ht = sum(dwn.*wn.*E,4)*dthd*dr;                                 % sum k^2*E over all directions                                                                               % spectral moment
            ht = sum(ht,3);                                                 % sum the prev sum over all frequencies to get the zeroth order moment (aka variance of sea surface)                                                                                                                               
            ht = 4*sqrt(abs(ht));                                           % signifigant wave height (from zeroth order moment of surface)                                                                       % sig wave height = 4*sqrt(spectral moment) for narrow-banded wave spectrum
            [sigH(iii,t),~] = max(max(ht));                                                                                                                     % return the largest signifigant wave height along the grid
            
            %m1 = sum((dwn.^2).*wn.*E,4)*dthd*dr;                            % first order moment m1 = integral f^2*S over all frequencies and direction
            m1 = sum((dwn).*(wn.^2).*E,4)*dthd*dr; 
            m1 = sum(m1,3);
            T0{iii,t} = ht./m1;                                             % avg time between zeroth-moment crossing (T0 = m0/m1)
            


            if showplots && rem(tplot,10) == 0    
               % Plot Signifigant Height
               [xplot,yplot] = meshgrid(1:m,1:n);
               figure;
               surf(xplot',yplot',ht,'EdgeColor','k')
               myc = colorbar;
               myc.Label.String = 'Sig H [m]';
               title(['Sig Wave Height for u = ',num2str(UUvec(iii)),' m/s'])
               frame = getframe(gcf);
               im{idx} = frame2im(frame);
               idx = idx + 1;
            end
% -- mean slope ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
           % integrate spectrum to find mean slope
           ms = sum(dwn.*wn.^3.*E,4)*dthd*dr;                                                                                                               % slope = angle water surface makes with flat surface
           ms = sum(ms,3);
           ms = sqrt(ms);                                                                                                                                   % standard deviation of water surface = sqrt(variance of water surface)

           if showplots && rem(tplot,10) == 0 
               figure(202);hold on;subplot(326);plot([1:m],ht(:,lati),'.-',[1:m],ms(:,lati),'--r');
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

           if showplots && rem(tplot,10) == 0
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
      
       fraction_time_completed = t/Tsteps;
       ttime=toc;
       seconds_so_far = toc;
       hours_to_complete = Tsteps/t*ttime/3600;
       fprintf('u = %.2f: fraction time completed: %.2f\n',UU,fraction_time_completed);
       %fprintf('Seconds so far: %i\n',ttime);
       %fprintf('Hours remaining: %.2f\n',hours_to_complete);
    
      
      
       eval(['save Titan/Titan_',int2str(file),'_',int2str(t),' E ht freqs oa Cd Cdf Cds Sds Sds_wc Sin Snl Sdt Sbf ms'])
       
       
       if t > 10 && sigH(iii,t-1)/sigH(iii,t) < tolH
           disp('Waves have reached 99% of maturity.')
           break
       elseif t > 1
           fprintf('t_n-1/t_n: %.6f\n',sigH(iii,t-1)/sigH(iii,t));
       end

   end % end of loop for Tsteps
  
   [minval,minind] = min(abs(UUvec-UU));
   disp(['Finished Wind Speed ' num2str(UU) ' m/s (' num2str(minind) ' Of ' num2str(numel(UUvec)) ' )']);
  
end % end of wind speed loop.
%% ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Make gif of results:
idx_end = idx;
filename = 'SigH.gif';                                                     % Specify the output file name
for idx = 1:idx_end-1
   [A,map] = rgb2ind(im{idx},256);
   if idx == 1
       imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
   else
       imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
   end
end
diary off; % close logging function
end
% ==============================================================================================================================================================================================================================================================================================================================================================================================================================
% User-Defined Functions:   
function k = wavekgt( f, H, g, T, R, tol )
%
%    k = wavekgt( f, H, g, T, R, tol )
%
%    Returns a column vector of wavenumbers [k] (1/m) at frequencies [f] (Hz)
%    using the linear wave dispersion relation for water depth H (m).  tol
%    is an optional error bound on the predicted f's (the default is 1e-4).
%    T is surface tension in N/m = dynes/cm * .001: default is 0.074
%    R = water density in Kgm/m^3: default is
%    Note: [f] > 0 must increase monotonically from small to large.  This
%    routine is optimized for speed and works best when length(f) is large.
                            % some initial housekeeping
if nargin<6,  tol = 1e-4;  end         % default tolerance
% g = 9.80171;
c = 2*pi*sqrt(H/g);   
f = c * f(:);  Nf = length(f);         % non-dimensionalize f
B = T/(R*g*H*H);                       % non-dimensionalize T
k = zeros(Nf,1);
   k(1) = f(1)^2;                     % use deep water limit as initial guess
   for n=1:Nf
   dk = 1;
   if n>1
   k(n) = k(n-1);
   end
   while ( abs(dk) ) > tol         % Newton-Raphson iteration loop
       t = 1;
       if k(n) < 20
       t = tanh( k(n) ) ;
       end
       dk = -(f(n).^2 - k(n).*(1+B.*k(n).^2)*t)...
            ./ ( 3*B.*k(n).^2*t+t  + k(n)*(1+B.*k(n).^2)*( 1 - t.^2 ) ) ;
       k(n) = k(n) - dk ;
   end
   k(n) = abs( k(n) ) ;               % f(k) = f(-k), so k>0 == k<0 roots
end
N = min( find( k>50 ) ) ;              % inform user if err > tol
if isempty(N)  
    N = Nf;
end
err = abs( f(2:N) - sqrt( k(2:N).*tanh(k(2:N)) ) )./f(2:N) ;
%if max( err ) > tol,  fprintf('\n WAVEK: error exceeds %g \n', tol),  end
k = k / H ;                            % give k dimensions of 1/meter
end
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
   %
  
  
   len = length(U);
   z0 = 0.001;
   z1 = 0.002;
   % nu=0.156/10000;
   del = 0.00001;
  
   for j = 1:len
       m = 0;
       u = U(j);
       z1 = z0+2*del;
       %while abs(z0-z1) > del
      
       for k = 1:6
           m = m+1;
           z1 = z0;
           ust = 0.4*u/(log(z/z0));
           z0 = 0.132*nu/ust;
       end
       ustar(j) = ust;
   end
end



