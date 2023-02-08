clc
clear all
close all

% Calculates E(x,y,k,theta) from energy balance equation.

% The equation to solve is:
% E(x,y,k,th)_{n+1} = E(x,y,k,th)_n + delt * (-Cg cos(th)dE/dx -Cg sin(th)dE/dy + Sin + Snl - Sds).
% Dimensions of x, y, k, th are m, n, ,o, p.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL SET UP %%

%% global constants %%
kappa = 0.4; % Von-Karman constant

%% grid %%
m = 31; % number of gridpoints in x-direction
n = 15; % number of gridpoints in y-direction
p = 288; % number of angular (th) bins. Note: p must be factorable by 8 for octants.
o = 25; % numner pf frequency bins

%% delt:
mindelt = 0.0001; % Minimum delt to minimize oscillations of delt to very small values
maxdelt = 2000.0; % Limit of delt to prevent large values as wind ramps up at startup

%% frequency limits:
f1 = 0.05;
f2=35;

%% frequency bins
oa = 7; % top bin advected
ol = 1:oa; % bin numbers for long frequencies
os = oa+1:o; % bin numbers for short frequencies.
op = [p/2+1:p 1:p/2]; % opposite direction from p 
orm = p:-1:1; % reflected directions (100% reflection at y-boundaries)

% wind profile:
zref = 20; % reference height for wind speed [m]
z = 10; % height of measured wind speed [m]
wfac = 0.035; % winddrift fraction of Uz 

%% LAKE OF INTEREST %%
% Titan, Punga Mare
long = 29; % longitude
lati = 8; % latitude
nu = 0.0031/1e4; % kinematic viscotity [m2/s] for pure methane
sfcT = 0.018; % surface tension of lake [N/m]
rhow = 465; % density of liquid [kg/m3]

%% TITAN CONSTANTS %%
sfcpress = 1.5*[101300]; % surface pressure [Pa].
sfctemp = 92; % surface temperature [K]
kgmolwt = 0.028; % gm molecular weight [Kgm/mol]
RRR = 8.314; % Universal gas constant [J/K/mol]
rhoa = sfcpress*kgmolwt/(RRR*sfctemp); % air density [kg/m3]
g = 1.35; % gravity [m/s2]
nua = 0.0126/1e4; % atmospheric gas viscosity [m2/s]
rhorat = rhoa/rhow; % air/water density ratio

%% WIND CLIMATE %%
gust = 0; % Gust factor applied to wind at each time step
UUvec = [0.4:1:4.5]; % wind velocities of interest [m/s]
wind_direction = 0; % direction of incoming wind
Cd=1.2*ones(n,n); % coefficient of friction is uniform over entire grid at 1.2

%% wavenumber:
kutoff = 1000; % wavenumber cut-off due to Kelvin-Helmholtz
kcg = sqrt(g*rhow/sfcT); % wavenumber of slowest waves aka capillary-gravity wave minimum
kcgn = 1.0*kcg; % n power min is not shifted wrt kcg
kcga = 1.15*kcg; % a min is shifted above kcg

%% wave angle
waveang = ([0:p-1]-p/2+0.5)*(360/p)*(pi/180);
waveang = repmat(waveang,[o 1 m n]);
waveang=shiftdim(waveang,2);

% theta
th = ([0:p-1]-p/2+0.5)*(360/p)*(pi/180);
tth = ([0:p-1])*(360/p)*(pi/180); % angular difference between short waves and longer waves.

%% indices for refraction rotation.
cw = [p 1:p-1];
ccw = [2:p 1];

%% upwave indices for advective term.
xp = [1 1:m-1];
yp = [1 1:n-1];
xm = [2:m m];
ym = [2:n n];

%% not sure?
cp = find(cos(th) > 0);
cm = find(cos(th) < 0);
sp = find(sin(th) > 0);
sm = find(sin(th) < 0);

% reshape not sure
cth = repmat(cos(th),[o 1 m n]);
cth=shiftdim(cth,2);
sth = repmat(sin(th),[o 1 m n]);
sth=shiftdim(sth,2);
cth2 = repmat(cos(tth).^2,[o 1 m n]);
cth2=shiftdim(cth2,2);

%% facs (fraction?)
mss_fac = 240; % MSS adjustment to Sdiss. 2000 produces very satisfactory overshoot, 400 does not
Sds_fac = 1.0;% This is fac with variable power nnn in Sds.
Sdt_fac = 0.001;
Sbf_fac = 0.002;% Select for bottom roughness. 0.004 is for smooth sandy bottom as in GOM.

%% Geographic deltas
Delx = 1000.0; % Titan [m]
Dely = 1000.0; % Titan [m]
dely = Dely*ones(m,n,o,p);
delx = Delx*ones(m,n,o,p);
mindelx = min(squeeze(delx(1,:,1,1)));

%% Time to run model
numdays = 1;
Time = 10000; % full time to run
Tsteps =2; % number of times run?
slostart = 1; % what is this?

% for saving data
file = 1; % for filename 

% something
explim = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LAKE GEOMETRY

% LAKE EXAMPLE 1: Constant Depth
%D = 100*ones(m,n); % constant liquid depth [m]

% LAKE EXAMPLE 2: Conical Island (how is this a cone?)

ccm=36; % x-location of island
ccn=26; % y-location of island

D=zeros(m,n); % initialize depth array

for j1=1:m
    for j2=1:n
        D(j1,j2) = 2*abs((j1-ccm) + 1i*(j2-ccn)) - 9;
    end
end


D(D<=0) = 0; % all negative depths should be set to zero

% lateral boundaries are absorptive
D(:,1)=0;
D(1,:)= 0;
D(:,end)=0;
D(end,:)=0;

%% ocean currents
Uer=0*D + 0.0; % Eastward current [m/s]
Uei=0*D + 0.0; % Northward current [m/s]

%% PLOT OF BATHYMETRY
[xplot,yplot] = meshgrid(1:m,1:n);
figure;
surf(xplot',yplot',D,'EdgeColor','k')
c = colorbar;
c.Label.String = 'Liquid Depth [m]';
title('Lake Model Bathymetry')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Snl %%

% dominant frequency (I think)
dlnf=(log(f2)-log(f1))/(o-1);
f = exp(log(f1)+[0:o-1]*dlnf);
dom = 2*pi*dlnf.*f;

% diffusion in two frequencies (bf1, bf2) and in two directions (bt1, bt2)
%   bf1 + bf2 = 1
%   bt1 + bt2 = 1
bf1 = exp(-15.73*dlnf*dlnf);
bf2 = exp(-15.73*4*dlnf*dlnf); 
bf1a = bf1/(bf1 + bf2);
bf2 = bf2/(bf1 + bf2);
bf1 = bf1a;

% Compute Snl_fac (downshifting is not a function of A)
A = exp(dlnf);
B = A.*A;
fac = bf1.*(1-1/B) + bf2.*(1-1./B.^2);
Snl_fac = ((1.0942/A)^1.9757)/fac;

fprintf('Snl_fac: %.3f\n',Snl_fac);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sdiss %%

% initialize for computing wavenumbers at each depth in the x,y grid

wn(:,:,:) = ones(m,n,o);
nnn(:,:,:) = ones(m,n,o);
ann(:,:,:) = ones(m,n,o);

addpath 'C:\Users\schne\OneDrive\Desktop\Main\Work\MIT\Github_Repos\umwm_titan' % need to add wavegkt() function call

for jm = 1:m
    for jn = 1:n
        if D(jm,jn) > 0
            wn(jm,jn,:) = WAVEKgT(f,D(jm,jn),g,sfcT,rhow,1e-4); % wave number 
            nnn(jm,jn,:) = 1.2 + 1.3*(abs(2 - (1+3*(wn(jm,jn,:)./kcgn).^2)./(1+(wn(jm,jn,:)./kcgn).^2)).^2.0); % Power of Sds
            ann(jm,jn,:) = 0.04 + 41.96*(abs(2 - (1+3*(wn(jm,jn,:)./kcga).^2)./(1+(wn(jm,jn,:)./kcga).^2)).^4.0); % Power of Sds
        end
    end
end

% reshape the matrix
wn = repmat(wn,[1 1 1 p]);
nnn = repmat(nnn,[1 1 1 p]);
ann = repmat(ann,[1 1 1 p]);
nnn = (2.53/2.5).*nnn;
nnninv = 1./nnn;
D=repmat(D,[1 1 o p]);
f=repmat(f',[1 p m n]);
f=shiftdim(f,2);
dom=repmat(dom',[1 p m n]);
dom=shiftdim(dom,2);

% wave phase velocity
c=(2*pi*f)./wn;
c(D<=0) = 0; % set phase velocity for negative depth to zero

% wave celerity
Cg=zeros(size(c)); % initialize
Cg(D>0) = c(D>0)./2.*(1 + 2*wn(D>0).*D(D>0)./sinh(2*wn(D>0).*D(D>0)) + 2*sfcT.*wn(D>0)./rhow./(g./wn(D>0) + sfcT.*wn(D>0)./rhow)); % Group velocity for all waves [m/s]
Cg(D<=0) = 0;
cgmax = max(max(max(max(Cg))));

% wave frequency?
dwn = ones(m,n,o,p); 
dwn(D>0)=dom(D>0)./abs(Cg(D>0));

% half wavelength
l2=abs(c)./f/2; %wavelength/2;
lz=abs(c)./f/2/pi; % wavelength/2/pi: kz = 1 for drift current action;

% Set half wavelength terms > zref equal to zref 
l2(l2 > zref) = zref;
lz(lz > zref) = zref;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = zeros(m,n,o,p); % initializing the energy term in x,y,k,theta with zeros
% Set delt by Courant condition (to enforce stability during time evolution in the finite differencing scheme)
Ccg = Cg + Uer.*cth + Uei.*sth;
delt = 0.7 * mindelx / cgmax; % advection is limited by the Courant condition for stability


% save data to mat file
eval(['save Titan/Titan_',int2str(file),'_ref m n o p freqs f wn dwn D Uei Uer Cg c delx dely dthd waveang'])
file = file - 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iii = 1:numel(UUvec) % iterate through all wind velocities of interest
    
    file = file + 1; % used for saving data to a specific filename and number
    modt = 0; % something for the time iteration

    U = UUvec(ii).*ones(m,n); % initialize wind vector to have same velocity everywhere on the grid
    U_z = U + 0.005; % need this small addition to avoid division by zero

    windir =  repmat(wind_direction.*ones(m,n),[1 1 o p]); % direction vector map of incoming wind (here it is set to be all uniform direction)
    uj = find(U>11); % for some reason
    Cd(uj) = 0.49 +.065.*U(uj); % for some reason
    Cd=Cd./1000; % for some reason

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for t = 1:Tsteps

        sumt = 0; % used for iterating 
        tplot = -1; % what is this?
        Newdelt = []; % initialize something

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        while((Time - sumt) > 0)

            tplot = tplot + 1;

            % wind velocity profile above an interface
            if gust ~= 0 
                U = U.*(1 + gust*randn(size(U))); % adds a random gust to the wind vector field
            end
            Ul = real(U_z.*sqrt(Cd)); % wind friction velocity 
            Ustar = Ul; % save copy to use elsehwere
            U_10 = 2.5*Ustar.*log(10/z) + U_z; % wind velocity profile according to the Law of the Wall
            Ustarw = Ustar .* sqrt(rhorat); % something related to shear stress

            % reshape matrix
            Ul = repmat(Ul,[1 1 o p]);
            U = repmat(U_z,[1 1 o p]);
            Ul = 2.5*Ul.*log(l2/z) + U ;    

            %calculate wind input as a fraction of stress.
            %   Ul = (Ucos(th)-c-drift*cos(phi)).^
            Ul(D>0)=abs(Ul(D>0).*cos(windir(D>0)-waveang(D>0))-c(D>0)-Ud(D>0).*cos(windir(D>0)-waveang(D>0))-Uer(D>0).*cth(D>0)-Uei(D>0).*sth(D>0))...
                .*(Ul(D>0).*cos(windir(D>0)-waveang(D>0))-c(D>0)-Ud(D>0).*cos(windir(D>0)-waveang(D>0))-Uer(D>0).*cth(D>0)-Uei(D>0).*sth(D>0));

            % adjust input to lower value when waves overrun the wind (as wave speed approaches wind speed the energy extraction becomes less efficient)
            Ula = Ul; % saving a copy
            Ul = 0.11*U1; % for some reason
            Ul(Ul<0) = Ula(Ul<0)*0.03;% where is this fraction from? UGS
            Ul(cos(windir-waveang)<0) = Ula(cos(windir-waveang)<0)*0.015;% Waves against wind.

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Sin -- INPUT TERM

            % input to energy spectrum
            Sin = zeros(size(Ul)); % initialize input
            Sin(D>0) = rhorat*(Ul(D>0).*2*pi.*f(D>0).*wn(D>0)./(g+sfcT.*wn(D>0).*wn(D>0)./rhow));
            % limit energy going into spectrum once waves outrun the wind by using a heavyside function
            Heavyside = ones(size(E)); % ones everywhere until waves outrun the wind
            Heavyside(Sin < 0 & E < 1e-320) = 0;  % zeroed out           
            Sin = Heavyside.*Sin;
            % set Sin = 0 at the land boundaries 
            Sin(D<=0) = 0;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Sds -- DISSIPATION TERM
            %   (1) turbulent dissipation: Sdt = 4 nu_t k^2
            %   (2) bottom friction: Sbf
            for tj = 1:p
                short(:,:,:,tj) = sum(E.*cth2(:,:,:,rem([1:p]-tj+p,p)+1),4)*(360/p)*(pi/180); % something
            end
            short = (cumsum((wn.^3.*short.*dwn),3)-wn.^3.*short.*dwn); % integrated up

            % Calculate turbulent dissipation
            Sdt(:,:,:,:) = repmat(Ustar,[1 1 o p]);
            Sdt = Sdt_fac*sqrt(rhorat).*Sdt.*wn;
            % Calculate bottom friction
            Sbf = zeros(m,n,o,p); % initialize
            Sbf(D>0) = Sbf_fac*wn(D>0)./sinh(2*wn(D>0).*D(D>0));

            % Calculate full dissipation term
            Sds = zeros(size(Sin)); % initialize
            Sds(D>0) = abs(ann(D>0)*2*pi.*f(D>0).*(1+mss_fac*short(D>0)).^2.*(wn(D>0).^4.*E(D>0)).^(nnn(D>0)));
            % set Sds = 0 at land boundaries
            Sds(D<=0) = 0; % LINE 469

        end
    end

end
