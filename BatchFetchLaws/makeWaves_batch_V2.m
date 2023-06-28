function [job] = makeWaves_batch_V2(planet,liquid,model,wind,uniflow,Etc)
%% ==========================================================================================================================================================================================================================================================================
%% ==========================================================================================================================================================================================================================================================================
% MAKEWAVES calculates E(x,y,k,theta) for wave field using an energy balance between wind-input and multiple dissipation terms (see Donelan et al. 2012 Modeling Waves and Wind Stress).
%
% The equation to solve is:
%   E_{n+1} = E_{n} + del*(-Cg*cos(theta)*dE/dx - Cg*sin(theta)*dE/dy] + Sin + Snl - Sds
%
%
% Dimensions of x,y,k,theta are m,n,o,p
%
%
%   Arguments:
%       planet
%           nua             : atmospheric gas viscosity [m2/s]
%           gravity         : gravitational acceleration [m/s2]
%           surface_temp    : surface temperature [K]
%           surface_press   : surface pressure [Pa]
%       liquid
%           rho_liquid      : liquid density [kg/m3]
%           nu_liquid       : liquid kinematic viscocity [m2/s]
%           surface_tension : surface tension of liquid [N/m]
%       model
%           m               : number of grid cells in x-dimension
%           n               : number of grid cells in y-dimension
%           o               : number of frequency bins
%           p               : number of angular (theta) bins. Must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
%           gridX           : size of grid cell in x-dimension [m]
%           gridY           : size of grid cell in y-dimension [m]
%           mindelt         : minimum time step to minimize oscillations of delt to very small values [s]
%           maxdelt         : maximum time step to prevent large values as wind ramps up at startup [s]
%           tolH            : minimum change in SigH to stop model at (otherwise runs to full Tsteps) 
%           bathy_map       : m x n array of depth [m] (+ values = subsurface, - values = subaerial)
%           time_step       : maximum size of time step [s]
%           num_time_steps  : length of model run in terms of # of time steps
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
%       sigH                : largest signifigant wave height over the entire spatial array [m]
%       htgrid              : signifigant wave height for each grid cell [m]
%       E_each              : wave energy spectrum (x,y) in space and in (frequency,direction) space
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
tic
%% ==========================================================================================================================================================================================================================================================================
%% ==========================================================================================================================================================================================================================================================================
assert(rem(model.p,8)==0,'Model input parameter p must be factorable by 8.')
% -- prepare log file for commands -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dfile=[datestr(now,'mmddyyyy_HHMMSS'),'_FetchLaws.txt'];
diary(dfile);
RAII.diary = onCleanup(@() diary('off'));                                  % auto-closes logging function on error
% -- create output directory for results -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TitanResults = sprintf('%s\\Titan%s', pwd);
if ~exist(TitanResults, 'dir')  && Etc.savedata                            % make output directory 'Titan' if doesn't already exist
    mkdir(TitanResults);
else
   oldmatfiles = fullfile(TitanResults,Etc.name,'*.mat');                          % empties output directory from previous runs
   oldmatloc = dir(oldmatfiles);
   for kk = 1:length(oldmatloc)
       basemat = oldmatloc(kk).name;
       fullmat = fullfile(TitanResults,basemat);
       fprintf(1,'Deleting previous .mat files %s\n',fullmat);
       delete(fullmat);
   end
end
disp('================================================================')
disp(['Directional Wave Spectrum -- last updated: ' dir('makeWaves.m').date])
disp(['Wind Speed(s) to Run: ' regexprep(mat2str(wind.speed),{'\[', '\]', '\s+'}, {'', '', ','}) ' m/s']);
disp('================================================================')
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
dthd = 360/(model.p);                                                      % step size for angular direction [degrees]
% Set up geographic deltas
dely = model.gridY*ones(model.m,model.n,model.o,model.p);
delx = model.gridX*ones(model.m,model.n,model.o,model.p);
% Wind:
gust = 0;                                                                  % Gust factor applied to wind at each time step
zref = 20;                                                                 % height of reference wind speed, normally at 20m [m]
z = 10;                                                                    % height of measured wind speed [m]
wfac = 0.035;                                                              % winddrift fraction of Uz ffor U10m
% Ideal gas law: PV = nRT
% Densities:
rhoa = planet.surface_press*kgmolwt/(RRR*planet.surface_temp);             % air density [kg/m3]
rhorat=rhoa/liquid.rho_liquid;                                             % air-water density ratio.
% Wavenumber limits:
kutoff = 1000;                                                             % wavenumber cut-off due to Kelvin-Hemholtz instabilities (Donelan 2012, pg. 3)
kcg = sqrt(planet.gravity*liquid.rho_liquid/liquid.surface_tension);       % wavenumber of slowest waves defined by the capillary-gravity waves, from Airy dispersion: omega^2 =gktanh{k)
% modified wavenumbers
kcgn = 1.0*kcg;                                                            % n power min is not shifted with respect to kcg
kcga = 1.15*kcg;                                                           % a min is shifted above kcg.
% frequency limits:
f1 = 0.05;                                                                 % minimum frequency
f2 = 35;                                                                   % maximum frequency 
% create frequency limits for spectrum
dlnf=(log(f2)-log(f1))/(model.o-1);                                        % frequency step size for log normal distribution
f = exp(log(f1)+(0:model.o-1)*dlnf);                                       % frequencies for spectrum
dom = 2*pi*dlnf.*f;                                                        % angular frequency (w = 2pi*f)
freqs = f;                                                                 % save a copy of frequencies
% Frequency bins:
oa = 7;                                                                    % top bin advected (make wavelengths [vector wn] at oa 1/20th of the grid spacing)
ol = 1:oa;                                                                 % bins for long frequencies
os = oa+1:model.o;                                                         % bins for short frequencies
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
% waveangle [radians]
waveang = ((0:model.p-1)-model.p/2+0.5)*dthd*dr;                           % wave angle [radians]
th = ((0:model.p-1)-model.p/2+0.5)*dthd*dr;                                % angle of wave propogation (phi in Donelan 2012)  (0.5 added to avoid division by zero in rotation term)
cth=cos(th);                                                               % cosine of angle in wave propogation direction 
sth=sin(th);                                                               % sine of angle in wave propogation direction
dth=dthd*dr;                                                               % 2pi/p (small angle for integration [radians])
% compute cos^2 for calculation of mss vs direction.
tth = (0:model.p-1)*dthd*dr;                                               % angular difference between short waves and longer waves
cth2 = cos(tth).^2;                                                        % cosine square of angular difference between short waves and longer waves
% indices for refraction rotation:
cw = ([model.p 1:model.p-1]);                                              % clockwise indices
ccw = ([2:model.p 1]);                                                     % counterclockwise indices
% upwave indices for advective term:
xp = [1 1:model.m-1];
yp = [1 1:model.n-1];
xm = [2:model.m model.m];
ym = [2:model.n model.n];
% Index being used for advection:
cp = find(cth > 0);
cm = find(cth < 0);
sp = find(sth > 0);
sm = find(sth < 0);
% reshape matrices
cth = repmat(cth,[model.o 1 model.m model.n]);
cth=shiftdim(cth,2);
sth = repmat(sth,[model.o 1 model.m model.n]);
sth=shiftdim(sth,2);
cth2 = repmat(cth2,[model.o 1 model.m model.n]);
cth2=shiftdim(cth2,2);
%% -- Lake Geometry ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
D = model.bathy_map;                                                       % depth of liquid [m]
D(D<=0) = 0;                                                               % limit land elevations to 0 to avoid dD/dx, dD/dy errors in refraction calculation                                                        
D(:,1) = 0; D(1,:) = 0; D(:,end) = 0; D(end,:) = 0;                        % Set array boundary depths to 0 (absorbtive boundaries)
% plot the bathymetry
[xplot,yplot] = meshgrid(1:model.m,1:model.n);
if Etc.showplots
    figure;
    surf(xplot,yplot,D','EdgeColor','k')
    myc = colorbar;
    myc.Label.String = 'Liquid Depth [m]';
    title('Lake Model Bathymetry')
end
%% -- Fractions of terms -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Fractions applied to terms from fitting with experiment/observations:
mss_fac = 240;                                                             % MSS adjustment to Sdiss. 2000 produces very satisfactory overshoot, 400 does not. (A4 in Donelan 2012; eqn. 6 Donelan 2001)
Sdt_fac = 0.001;                                                           % fraction of Sdt that goes into spectrum (A4 in eqn. 20, Donelan 2012; A4 = 0.01?)
Sbf_fac = 0.002;                                                           % fraction of Sbf that goes into spectrum (0.004 is for smooth sandy bottom) (Gf in eqn. 22, Donelan 2012) 
%% -- wavemumber and Power of Sds ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
wn(:,:,:) = ones(model.m,model.n,model.o);
nnn(:,:,:) = ones(model.m,model.n,model.o);
ann(:,:,:) = ones(model.m,model.n,model.o);
m = model.m; n = model.n; g = planet.gravity;
sfcT = liquid.surface_tension; rhow = liquid.rho_liquid; 
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
wn = repmat(wn,[1 1 1 model.p]);
nnn = repmat(nnn,[1 1 1 model.p]);
ann = repmat(ann,[1 1 1 model.p]);
% conversions for nnn?
nnn = (2.53/2.5).*nnn;
nnninv = 1./nnn;
% reshape the matrix 
D = repmat(D,[1 1 model.o model.p]);
f = repmat(f',[1 model.p model.m model.n]);
f = shiftdim(f,2);
dom = repmat(dom',[1 model.p model.m model.n]);
dom = shiftdim(dom,2);
%% -- ocean currents ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Uer = uniflow.East*D + 0.0;                                                                                                                            % Eastward current, m/s
Uei = uniflow.North*D + 0.0;                                                                                                                           % Northward current, m/s
%% -- wave speed and group velocity --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
c = (2*pi*f)./wn;                                                                                                                                                                                                    % phase speed
c(D<=0) = 0;                                                                                                                                                                                                         % phase speed on land is set to zero
Cg = zeros(size(c));                                                                                                                                                                                                 % initialize group speed
Cg(D>0) = c(D>0)./2.*(1 + 2*wn(D>0).*D(D>0)./sinh(2*wn(D>0).*D(D>0)) + 2*liquid.surface_tension.*wn(D>0)./liquid.rho_liquid./(planet.gravity./wn(D>0) + liquid.surface_tension.*wn(D>0)./liquid.rho_liquid));        % Group velocity for all waves (Kinsman "Wind Waves: Their Generation and Propogation on the Ocean Surface")
dwn = ones(model.m,model.n,model.o,model.p);                                                                                                                                                                         % initalize dominant wavenumber (c = dw/dk)
dwn(D>0) = dom(D>0)./abs(Cg(D>0));                                                                                                                                                                                   % remove any values on land for dominant wavenumber (c = dw/dk)
Cg(D<=0) = 0;                                                                                                                                                                                                        % zero all the group velocities on land
l2=abs(c)./f/2;                                                                                                                                                                                                      % wavelength/2
lz=abs(c)./f/2/pi;                                                                                                                                                                                                   % wavelength/2/pi: kz = 1 for drift current action
% Set l2, lz > zref equal to zref   
l2(l2 > zref) = zref;
lz(lz > zref) = zref;
E = zeros(model.m,model.n,model.o,model.p);                                                                                                                                                                          % initializing the energy term in x,y,k,theta with zeros
cgmax = max(max(max(max(Cg))));                                                                                                                                                                                      % fastest group velocity
mindelx = min(squeeze(delx(1,:,1,1)));                                                                                                                                                                               % smallest spatial grid spacing for domain of influence
waveang = repmat(waveang,[model.o 1 model.m model.n]);
waveang=shiftdim(waveang,2);
if Etc.savedata
    m = model.m; n = model.n; o = model.o; p = model.p;
    save([TitanResults,'/New_Reference_',Etc.name],'m','n','o','p','freqs','f','wn','dwn','D','Uei','Uer','Cg','c','delx', 'dely','dthd','waveang','dr')
end
file = file - 1;
%% -- loop through wind speeds --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
wind.speed = [0 1 2 3 4];
idx = 1;                                                                   % for frame for gif
sigH = zeros(numel(wind.speed),length(1:model.num_time_steps));            % initialize sigH for returning 
htgrid = cell(1,numel(wind.speed));                                        % initialize htgrid for returning
E_each =  cell(numel(wind.speed),length(1:model.num_time_steps));          % initialize E-spectrum for returning
%% -- batch processing wind speeds ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
range = 1:size(wind.speed,1);
numwinds = size(range,1);
job{numwinds} = [];

for Wind = range
    job{Wind} = batch('winds_batch',3,{planet,liquid,model,wind,Etc,ann,bf1,bf2,c,ccw,Cg,cgmax,cm,cp,cth,cth2,cw,D,delx,dely,dr,dth,dthd,dwn,E,E_each,f,file,freqs,gust,htgrid,i,idx,kappa,kutoff,l2,lz,mindelx,mss_fac,nnn,nnninv,oa,ol,os,rhoa,rhorat,Sbf_fac,Sdt_fac,sigH,sm,Snl_fac,sp,sth,TitanResults,waveang,wfac,wn,xm,xp,ym,yp,z});% [m]
            pause(.01)
end
