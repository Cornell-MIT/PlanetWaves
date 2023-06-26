function [job] = slope_effect_batch
% make plots of signifigant wave height for slopes

%% LOAD DATA:
% top row = methane alkane fraction (molar fraction of the alkanes (CH4+C2H6) that is methane)
% first column = temperature [K]
format long
density = load('../compositions/density.csv'); % [kg/m3]
kinematic_visc = load('../compositions/kinematic_visc.csv'); % [cm2/s]
dynamic_visc = load('../compositions/dynamic_visc.csv'); % [Pa*s]
surf_ten = load('../compositions/surface_tension.csv'); % [N/m]
methane_frac = load('../compositions/CH4_molefraction.csv'); % methane
ethane_frac = load('../compositions/C2H6_molefraction.csv'); % ethane
nitrogen_frac = load('../compositions/N2_molefraction.csv'); % nitrogen

%% INPUT PARAMETERS:
% (1) PLANET CONDITIONS
%   (a) TITAN
Titan.nua = 0.0126/1e4;                                                    % Titan atmospheric gas viscocity [m2/s]
Titan.gravity = 1.352;                                                     % Titan Gravity [m/s2]
Titan.surface_temp = 92;                                                   % Titan Surface Temperature [K]
Titan.surface_press = 1.5*101300;                                          % Titan Surface Pressure [Pa]
Titan.liquid.Hydrocarbon.rho_liquid = 540;                                 % Hydrocarbon liquid density [kg/m3]
Titan.liquid.Hydrocarbon.nu_liquid = 3e-7;                                 % Hydrocarbon liquid Viscocity [m2/s]
Titan.liquid.Hydrocarbon.surface_tension = 0.018;                              % Hydrocarbon liquid Surface Tension [N/m]

%   (a.1) VARYING METHANE:ETHANE:NITROGEN COMPOSITIONS @92 K               % Source: steckloff et al., 2020
temp_location = find(methane_frac(:,1)==Titan.surface_temp);               % Find row with temperature
methane = 0:0.1:1;                                                         % Fraction of methane
ethane = 1:-0.1:0;                                                         % Fraction of ethane
for i = 1:length(methane)
    nitrogen(i)=nitrogen_frac(temp_location,...                            % Fraction of nitrogen
        max(find(abs(nitrogen_frac(1,:) - methane(i))<0.001)));
    name(i) = string(sprintf('m%03.0fe%03.0fn%03.0f',...                   % Creating structural array with names
        100*methane(i),100*ethane(i),100*nitrogen(i)));
    Titan.liquid.(name(i)).rho_liquid = density(temp_location,...          % Liquid density [kg/m3]
        max(find(abs(density(1,:) - methane(i))<0.001))); 
    Titan.liquid.(name(i)).nu_liquid = kinematic_visc(temp_location,...    % Liquid viscocity [m2/s]
        max(find(abs(kinematic_visc(1,:) - methane(i))<0.001)))/10000;      
     Titan.liquid.(name(i)).surface_tension = surf_ten(temp_location,...   % Liquid surface tension [N/m]
        max(find(abs(surf_ten(1,:) - methane(i))<0.001))); 
end
                                                      

%   (a.2) VARYING LAKE COMPOSITIONS @92 K
%   (a.2.1) TITAN: ONTARIO LACUS
%   51:38:11 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
Ontario_methane = 0.57;                                                    % Fraction of methane
Ontario_ethane = 0.43 ;                                                    % Fraction of ethane
Ontario_nitrogen = nitrogen_frac(temp_location,...                         % Fraction of nitrogen
        max(find(abs(nitrogen_frac(1,:) - Ontario_methane)<0.001)));
Titan.liquid.Ontario.rho_liquid = density(temp_location,...                % Ontario Lacus liquid density [kg/m3]
        max(find(abs(density(1,:) - Ontario_methane)<0.001))); 
Titan.liquid.Ontario.nu_liquid = kinematic_visc(temp_location,...          % Ontario Lacus liquid Viscocity [m2/s]
        max(find(abs(kinematic_visc(1,:) - Ontario_methane)<0.001)))/10000;      
Titan.liquid.Ontario.surface_tension = surf_ten(temp_location,...          % Ontario Lacus liquid Surface Tension [N/m]
        max(find(abs(surf_ten(1,:) - Ontario_methane)<0.001)));                                                       

%   (a.2.2) TITAN: PUNGA MARE
%   80:0:20 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
Punga_methane = 1;                                                         % Fraction of methane
Punga_ethane = 0 ;                                                         % Fraction of ethane
Punga_nitrogen = nitrogen_frac(temp_location,...                           % Fraction of nitrogen
        max(find(abs(nitrogen_frac(1,:) - Punga_methane)<0.001)));
Titan.liquid.Punga.rho_liquid = density(temp_location,...                  % Punga Mare liquid density [kg/m3]
        max(find(abs(density(1,:) - Punga_methane)<0.001))); 
Titan.liquid.Punga.nu_liquid = kinematic_visc(temp_location,...            % Punga Mare liquid Viscocity [m2/s]
        max(find(abs(kinematic_visc(1,:) - Punga_methane)<0.001)))/10000;      
Titan.liquid.Punga.surface_tension = surf_ten(temp_location,...            % Punga Mare liquid Surface Tension [N/m]
        max(find(abs(surf_ten(1,:) - Punga_methane)<0.001)));                                                           

%   (a.2.3) TITAN: LIGIA MARE (MAX METHANE)
%   100:0:0 percent methane:ethane:nitrogen from LeGall 2016
LigiaMax_methane = 1;                                                      % Fraction of methane
LigiaMax_ethane = 0 ;                                                      % Fraction of ethane
LigiaMax_nitrogen = nitrogen_frac(temp_location,...                        % Fraction of nitrogen
        max(find(abs(nitrogen_frac(1,:) - LigiaMax_methane)<0.001)));
Titan.liquid.LigiaMax.rho_liquid = density(temp_location,...               % Ligia Mare (max) liquid density [kg/m3]
        max(find(abs(density(1,:) - LigiaMax_methane)<0.001))); 
Titan.liquid.LigiaMax.nu_liquid = kinematic_visc(temp_location,...         % Ligia Mare (max) liquid Viscocity [m2/s]
        max(find(abs(kinematic_visc(1,:) - LigiaMax_methane)<0.001)))/10000;      
Titan.liquid.LigiaMax.surface_tension = surf_ten(temp_location,...         % Ligia Mare (max) liquid Surface Tension [N/m]
        max(find(abs(surf_ten(1,:) - LigiaMax_methane)<0.001)));                                                            

%   (a.2.4) TITAN: LIGIA MARE (MIN METHANE)
%   50:39:10.73 percent methane:ethane:nitrogen from LeGall 2016
LigiaMin_methane = 0.56;                                                   % Fraction of methane
LigiaMin_ethane = 0.44 ;                                                   % Fraction of ethane
LigiaMin_nitrogen = nitrogen_frac(temp_location,...                        % Fraction of nitrogen
        max(find(abs(nitrogen_frac(1,:) - LigiaMin_methane)<0.001)));
Titan.liquid.LigiaMin.rho_liquid = density(temp_location,...               % Ligia Mare (min) liquid density [kg/m3]
        max(find(abs(density(1,:) - LigiaMin_methane)<0.001))); 
Titan.liquid.LigiaMin.nu_liquid = kinematic_visc(temp_location,...         % Ligia Mare (min) liquid Viscocity [m2/s]
        max(find(abs(kinematic_visc(1,:) - LigiaMin_methane)<0.001)))/10000;      
Titan.liquid.LigiaMin.surface_tension = surf_ten(temp_location,...         % Ligia Mare (min) liquid Surface Tension [N/m]
        max(find(abs(surf_ten(1,:) - LigiaMin_methane)<0.001)));                                                            
                                                  

%   (a.2.5) TITAN: KRAKEN MARE
%   70:16:14 percent methane:ethane:nitrogen from Poggiali et al., 2020
Kraken_methane = 0.81;                                                     % Fraction of methane
Kraken_ethane = 0.19 ;                                                     % Fraction of ethane
Kraken_nitrogen = nitrogen_frac(temp_location,...                          % Fraction of nitrogen
        max(find(abs(nitrogen_frac(1,:) - Kraken_methane)<0.001)));
Titan.liquid.Kraken.rho_liquid = density(temp_location,...                 % Kraken liquid density [kg/m3]
        max(find(abs(density(1,:) - Kraken_methane)<0.001))); 
Titan.liquid.Kraken.nu_liquid = kinematic_visc(temp_location,...           % Kraken liquid Viscocity [m2/s]
        max(find(abs(kinematic_visc(1,:) - Kraken_methane)<0.001)))/10000;      
Titan.liquid.Kraken.surface_tension = surf_ten(temp_location,...           % Kraken liquid Surface Tension [N/m]
        max(find(abs(surf_ten(1,:) - Kraken_methane)<0.001)));                                             

%   (b) EARTH
Earth.liquid.Water.rho_liquid = 997;                                       % Water liqid density [kg/m3]         
Earth.liquid.Water.nu_liquid = 1e-6;                                       % Water liquid viscocity [m2/s]
Earth.liquid.Water.surface_tension = 0.072;                                    % Water liquid surface ension [N/m]
Earth.gravity = 9.81;                                                      % Earth Gravity [m/s2]
Earth.surface_temp = 273;                                                  % Earth Surface Temperature [K]
Earth.surface_press = 1*101300;                                            % Earth Surface Pressure [Pa]
Earth.nua = 1.338/1e5;                                                     % Earth atmospheric gas viscocity [m2/s]

%% (2) MODEL GEOMETRY
Model.m = 30;                                                              % Number of Grid Cells in X-Dimension
Model.n = 30;                                                              % Number of Grid Cells in Y-Dimension
Model.o = 25;                                                              % Number of Frequency bins
Model.p = 288;                                                             % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.long = 15;                                                           % longitude grid point for sampling during plotting 
Model.lat = 15;                                                            % latitude grid point for sampling during plotting
Model.gridX = 1000.0;                                                      % Grid size in X-dimension [m]
Model.gridY = 1000.0;                                                      % Grid size in Y-dimension [m]
Model.mindelt = 0.0001;                                                    % minimum time step
Model.maxdelt = 2000.0;                                                    % maximum time step
Model.time_step = 100;                                                    % MAXIMUM SIZE OF TIME STEP [S]
Model.num_time_steps = 100;                                                % LENGTH OF MODEL RUN (IN TERMS OF # OF TIME STEPS)
Model.tolH = NaN;                                                          % TOLERANCE THRESHOLD FOR MATURITY 

% % define a bathymetry with a constant slope
% bathy_map = ones(Model.m,Model.n);
% for i = 1:Model.m
%     bathy_map(i,:) = 50*(i/Model.m);
% end
% Model.bathy_map = bathy_map;                                             % Bathymetry of model basin [m]


% define a bathymetery map that is constant and deep
Model.bathy_map = 100.*ones(Model.m,Model.n);                              % Bathymetry of model basin [m]


%% (3) NEAR-SURFACE WIND CONDITIONS
Wind.dir = 0;                                                              % direction of incoming wind [radians]
%Wind.speed = log10(logspace(0,5,50));                                      % wind speed [m/s]
Wind.speed = [0 1 2 3 4];
%% (4) Unidirectional currents
Uniflow.East = 0;                                                          % eastward unidirectional current [m/s]
Uniflow.North = 0;                                                         % northward unidirectional current [m/s]

%% (5) HOUSEKEEPING
Composition = string(fieldnames(Titan.liquid));
range = 1:size(Composition,1);
numcomps = size(range,1);
job{numcomps} = [];

Etc.showplots = 0;
Etc.savedata = 0;
Etc.showlog = 0;
%% (6) RUN MODEL

% Modeling Titan Waves
    for liq = range
        Etc.name = convertStringsToChars(Composition(liq));
        %makeWaves_batch(Titan,Titan.liquid.(Composition(liq)),Model,Wind,Uniflow,Etc);   
    
        job{liq} = batch('makeWaves',3,{Titan,Titan.liquid.(Composition(liq)),Model,Wind,Uniflow,Etc});% [m]
            pause(.01)
    
    end
end

%S=load("mydata.mat");
%[data,hdr]=deal(S.data,S.hdr);




