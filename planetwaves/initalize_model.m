function [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_name,time_to_run,wind_dir,depth_profile,buoy_loc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function will populate model parameters with default suggested values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
%       planet_name     : name of planet [string] e.g. = 'Titan'
%       time_to_run     : number of time steps [#]
%       wind_dir        : direction of wind [radians]
%       depth_profile   : 2D array of depths [m]
%       buoy_loc        : grid location (x,y) of buoy for measurement
%   OUTPUT:
%       Planet          : Class object containing planetary conditions
%       Model           : Class object containing model geometry 
%       Wind            : Class object containing wind strength and direction
%       Uniflow         : Class object containing unidirectional North and Eastward flow within basin
%       Etc             : Class object to plot intermediary results and save run data in .mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ATM_2_PASCAL = 101325;                                                  % Conversion bw atms and pascals


Planet.name                        = planet_name;                       % string name of planet
size_lake                          = size(depth_profile);               % size of array                      
Model.LonDim                       = size_lake(2);                      % Number of Grid Cells in X-Dimension (col count)
Model.LatDim                       = size_lake(1);                      % Number of Grid Cells in Y-Dimension (row count)

if ( (buoy_loc(1) == 1 || buoy_loc(2) == 1 || buoy_loc(1) == Model.LonDim || buoy_loc(2) == Model.LatDim))
    error('Buoy cannot be located at the edge of the grid. Choose different buoy_loc')
end

% -----------------------------------------------------------------------------------------------------------------------
% Array sizes and sampling
Model.Fdim                         = 50;                                % Number of Frequency bins
Model.Dirdim                       = 72;                                % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.long                         = buoy_loc(1);                       % longitude grid point for sampling during plotting
Model.lat                          = buoy_loc(2);                       % latitude grid point for sampling during plotting
% -----------------------------------------------------------------------------------------------------------------------
% ANENOMETER HEIGHT AND LIQUID DEPTH
Model.z_data                       = 10;                                % elevation of wind measurement [m]
Model.bathy_map                    = depth_profile;                     % bathymetric map [m]
% -----------------------------------------------------------------------------------------------------------------------
% Frequency Arrays
Model.min_freq                     = 0.05;                              % minimum frequency to model
Model.max_freq                     = 35;                                % maximum frequency to model
% -----------------------------------------------------------------------------------------------------------------------
% cut offs for diagnostic portion or maturity of field
%Model.cutoff_freq = round((20/35)*Model.Fdim);                         % cutoff frequency bin from diagnostic to advection -- if set too low can lead to numerical ringing
Model.tolH                        = NaN;                                % tolerance threshold for maturity
% -----------------------------------------------------------------------------------------------------------------------
% WINDS
Wind.dir                          = wind_dir;                           % wind direction (from -> to) [radians]
% -----------------------------------------------------------------------------------------------------------------------
% CURRENTS
Uniflow.East                      = 0;                                  % eastward unidirectional current [m/s]
Uniflow.North                     = 0;                                  % northward unidirectional current [m/s]
% -----------------------------------------------------------------------------------------------------------------------
% show plots during runs and saving data (if yes, 1)
Etc.showplots                     = 0;
Etc.savedata                      = 0;
% -----------------------------------------------------------------------------------------------------------------------
% Empirical-ish values for wave equations
Model.tune_A1                     = 0.11;                               % wind sea (eq. 12 Donelan+2012)
Model.tune_mss_fac                = 360;
Model.tune_Sdt_fac                = 0.001;
Model.tune_Sbf_fac                = 0.002;
Model.tune_cotharg                = 0.2;
Model.tune_n                      = 2.4;
Model.explim                      = 0.1;
% -----------------------------------------------------------------------------------------------------------------------
% MODEL TIME
% sub-time step parameters for time evolution 
Model.mindelt                     = 0.0001;                             % seconds
Model.maxdelt                     = 2000.0;                             % seconds
Model.time_step                   = 60;                                 % time step size [seconds]
Model.num_time_steps              = time_to_run;                        % total time [seconds]
% -----------------------------------------------------------------------------------------------------------------------
% molecular weights
N2_atm_kgmolwt                    = 0.028;
CO2_atm_kgmolwt                   = 0.044;
% experimental values from chemistry database: dippr.aiche.org DATABASE
% (1) Liquid properties:
%       WATER
liquid_H2O_density_288K           = 998.21;                             % [kg/m3] at 288K (DIPPR)
liquid_H2O_dyn_viscosity_288K     = 0.0011539;                          % [Pa.s] at 288 K (DIPPR)
liquid_H2O_surface_tension_288K   = 0.073836;                           % [N/m] at 288 K (DIPPR)
%       LIQUID NITROGEN
liquid_N2_density_77p5K           = 807.22;                             % [kg/m3] at 77.5K (DIPPR)  
liquid_N2_dyn_viscosity_77p5K     = 0.00015603;                         % [Pa.s] at 77.5K (DIPPR) 
liquid_N2_surface_tension_77p5K   = 0.0088249;                          % [N/m] at 77.5K (DIPPR)

% (2) Gas properties:
%       NITROGEN ATMOSPHERE
vapor_N2_dyn_viscosity_288K       = 0.000017253;                        % [Pa.s] at 288 K (DIPPR)  (e.g., Earth)
vapor_N2_dyn_viscosity_77p5K      = 0.0000054170;                       % [Pa.s] at 77.5 K (DIPPR) (e.g., past Titan)
vapor_N2_dyn_viscosity_92K        = 0.0000064320;                       % [Pa.s] at 92 K (DIPPR) (e.g., Titan)
%       CO2 ATMOSPHERE
vapor_CO2_dyn_viscosity_288K      = 0.000014482;                        % [Pa.s] at 288 K (DIPPR) (e.g., past Mars, like in Hort & Weitz 2001, which used 1.2e-5 Pa.s)
vapor_CO2_dyn_viscosity_1500K     = 0.000052032;                        % [Pa.s] at 1500 K (DIPPR)  
% -----------------------------------------------------------------------------------------------------------------------

if strcmp(planet_name,'Earth')

    % EARTH CONDITIONS
    Planet.rho_liquid            = liquid_H2O_density_288K;                            % Liquid density [kg/m3]                               
    Planet.nu_liquid             = liquid_H2O_dyn_viscosity_288K/Planet.rho_liquid;    % Liquid (kinematic) Viscosity (kinematic viscocity = dynamic viscocity / density of liquid) [m2/s]                        
    Planet.gravity               = 9.81;                                               % Gravity [m/s2]
    Planet.surface_temp          = 288;                                                % Surface Temperature [K]
    Planet.surface_press         = 1*ATM_2_PASCAL;                                     % Surface Pressure [Pa]      
    Planet.surface_tension       = liquid_H2O_surface_tension_288K;                    % Liquid Surface Tension [N/m] (DIPPR)
    Planet.kgmolwt               = N2_atm_kgmolwt;                                     % gram molecular weight [Kgm/mol] (N2)  
    Planet.rhoa                  = ideal_gas_law(Planet);                              % Atm density from ideal gas [kg/m3]
    Planet.nua                   = vapor_N2_dyn_viscosity_288K/Planet.rhoa;            % Atm gas (kinematic) viscosity (kinematic viscosity = dynamic viscosity / density of air) [m2/s] 
    
elseif strcmp(planet_name,'Mars-low')

    % MARS CONDITIONS AT JEZERO
    Planet.rho_liquid            = liquid_H2O_density_288K;                            % liquid water density at 288 K       
    Planet.nu_liquid             = liquid_H2O_dyn_viscosity_288K/Planet.rho_liquid;    % liquid (kinematic) Viscosity [m2/s]
    Planet.gravity               = 3.71;                                               % Mars gravity
    Planet.surface_temp          = 288;                                                % Earth-like warmth for liquid water to exist
    Planet.surface_press         = 50000;                                              % low pressure case
    Planet.surface_tension       = liquid_H2O_surface_tension_288K;                    % liquid water surface tension at 288 K (DIPPR)
    Planet.kgmolwt               = CO2_atm_kgmolwt;                                    % (CO2)
    Planet.rhoa                  = ideal_gas_law(Planet);                              % Atm density from ideal gas [kg/m3]
    Planet.nua                   = vapor_CO2_dyn_viscosity_288K/Planet.rhoa;           % CO2 vapor (kinematic) viscosity [m2/s]        

elseif strcmp(planet_name,'Mars-high')

    % MARS CONDITIONS AT JEZERO
    Planet.rho_liquid            = liquid_H2O_density_288K;                            % liquid water density at 288 K (DIPPR)         
    Planet.nu_liquid             = liquid_H2O_dyn_viscosity_288K/Planet.rho_liquid;    % liquid water viscocity at 288 K (DIPPR) 
    Planet.gravity               = 3.71;                                               % Mars gravity
    Planet.surface_temp          = 288;                                                % Earth-like warmth for liquid water to exist
    Planet.surface_press         = 4*50000;                                            % high pressure case
    Planet.surface_tension       = liquid_H2O_surface_tension_288K;                    % liquid water surface tension at 288 K (DIPPR)
    Planet.kgmolwt               = CO2_atm_kgmolwt;                                    % (CO2)
    Planet.rhoa                  = ideal_gas_law(Planet);                              % Atm density from ideal gas [kg/m3] 
    Planet.nua                   = vapor_CO2_dyn_viscosity_288K/Planet.rhoa;           % CO2 vapor (kinematic) viscosity [m2/s]  

elseif strcmp(planet_name,'Titan-N2')
    % TITAN CONDITIONS IN PAST/FUTURE WITH NO METHANE TO KEEP MOON WARM
    % liquid nitrogen (based on methane depletion case from Charnay+2014) for surface albedo of 0.3 with non-radiative clouds (Table 2)
    Planet.rho_liquid           = liquid_N2_density_77p5K;                            % liquid N2 density at 77.5 K (DIPPR)   
    Planet.nu_liquid            = liquid_N2_dyn_viscosity_77p5K/Planet.rho_liquid;    % liquid N2 viscocity at 77.5 (DIPPR)  
    Planet.gravity              = 1.352;                                              % Titan gravity
    Planet.surface_temp         = 77.5;                                               % Charney+2014                           
    Planet.surface_press        = 0.89*ATM_2_PASCAL;                                  % Charney+2014
    Planet.surface_tension      = liquid_N2_surface_tension_77p5K;                    % liquid N2 surface tension at 77.5 K (DIPPR)   
    Planet.kgmolwt              = N2_atm_kgmolwt;                                     % (N2)                                         
    Planet.rhoa                 = ideal_gas_law(Planet);                              % Atm density from ideal gas [kg/m3]
    Planet.nua                  = vapor_N2_dyn_viscosity_77p5K/Planet.rhoa;           % N2 vapor (kinematic) viscosity [m2/s]

elseif strcmp(planet_name,'Titan-OntarioLacus')
    % TITAN CONDITIONS at Ontario Lacus (ethane-rich)
    % 47% CH4, 40% C2H6, and 13% N2 [Mastrogiuseppe+2018] (values taken from Steckloff TitanPool for Methane Alkaline Fraction 0.47 @ 92K)
    Planet.rho_liquid          = 588.15;                                              % TITANPOOL   
    Planet.nu_liquid           = 8.084e-7;                                            % TITANPOOL   
    Planet.gravity             = 1.352;                                               % Titan gravity
    Planet.surface_temp        = 92;                                                  % ~avg Titan surface temp
    Planet.surface_press       = 1.5*ATM_2_PASCAL;                                    % ~avg Titan surface pressure
    Planet.surface_tension     = 0.032766;                                            % TITANPOOL
    Planet.kgmolwt             = N2_atm_kgmolwt;                                      % (N2)
    Planet.rhoa                = ideal_gas_law(Planet);                               % Atm density from ideal gas [kg/m3] (compare to 0.005 g/cm3 in Hayes+2013, Table 1) 
    Planet.nua                 = vapor_N2_dyn_viscosity_92K/Planet.rhoa;              % N2 vapor (kinematic) viscosity [m2/s]

elseif strcmp(planet_name,'Titan-LigeiaMare')
    % TITAN CONDITIONS at Ligea Mare (methane-rich)
    % 71% CH4, 12% C2H6, and 17% N2 [Mastrogiuseppe+2016] (values taken from Steckloff TitanPool for Methane Alkaline Fraction 0.71 @ 92K)
    Planet.rho_liquid         = 551.06;                                               % TITANPOOL  
    Planet.nu_liquid          = 5.685e-7;                                             % TITANPOOL   
    Planet.gravity            = 1.352;                                                % Titan gravity
    Planet.surface_temp       = 92;                                                   % ~avg Titan surface temp
    Planet.surface_press      = 1.5*ATM_2_PASCAL;                                     % ~avg Titan surface pressure 
    Planet.surface_tension    = 0.028363;                                             % TITANPOOL
    Planet.kgmolwt            = N2_atm_kgmolwt;                                       % (N2)
    Planet.rhoa               = ideal_gas_law(Planet);                                % Atm density from ideal gas [kg/m3]
    Planet.nua                = vapor_N2_dyn_viscosity_92K/Planet.rhoa;               % N2 vapor (kinematic) viscosity [m2/s]
                                  
elseif strcmp(planet_name,'Titan-CH4N2')
    % TITAN CONDITIONS: Only Methane and Nitrogen
    % [Mastrogiuseppe+2016] (values taken from Steckloff TitanPool for Methane Alkaline Fraction 1.0 @ 92K)
    Planet.rho_liquid        = 510.65;                                                % TITANPOOL   
    Planet.nu_liquid         = 3.1827e-7;                                             % TITANPOOL   
    Planet.gravity           = 1.352;                                                 % Titan gravity
    Planet.surface_temp      = 92;                                                    % ~avg Titan surface temp
    Planet.surface_press     = 1.5*ATM_2_PASCAL;                                      % ~avg Titan surface pressure
    Planet.surface_tension   = 0.016606;                                              % TITANPOOL
    Planet.kgmolwt           = N2_atm_kgmolwt;                                        % (N2)
    Planet.rhoa              = ideal_gas_law(Planet);                                 % Atm density from ideal gas [kg/m3]
    Planet.nua               = vapor_N2_dyn_viscosity_92K/Planet.rhoa;                % N2 vapor (kinematic) viscosity [m2/s]
                            
elseif strcmp(planet_name,'Titan-CH3H8N2')
    % TITAN CONDITIONS: Only Ethane and Nitrogen
    % [Mastrogiuseppe+2016] (values taken from Steckloff TitanPool for Methane Alkaline Fraction 0.0 @ 92K)
    Planet.rho_liquid        = 654.69;                                                % TITANPOOL   
    Planet.nu_liquid         = 1.6501e-6;                                             % TITANPOOL    
    Planet.gravity           = 1.352;                                                 % Titan gravity
    Planet.surface_temp      = 92;                                                    % ~avg Titan surface temp
    Planet.surface_press     = 1.5*ATM_2_PASCAL;                                      % ~avg Titan surface pressure
    Planet.surface_tension   = 0.03317;                                               % TITANPOOL
    Planet.kgmolwt           = N2_atm_kgmolwt;                                        % (N2)
    Planet.rhoa              = ideal_gas_law(Planet);                                 % Atm density from ideal gas [kg/m3]
    Planet.nua               = vapor_N2_dyn_viscosity_92K/Planet.rhoa;                % N2 vapor (kinematic) viscosity [m2/s]

elseif strcmp(planet_name,'Kepler-1649-b') % Exo-Venus
    % Sulfuric Acid
    Planet.rho_liquid        = 1838.1;                                                % liquid sulfuric acid density at 288 K (DIPPR)                                                                         
    Planet.nu_liquid         = 3.4e-5;                                                % liquid sulfuric acid viscocity at 288 K (DIPPR)
    Planet.gravity           = 8.87;                                                  % Venus-like gravity
    Planet.surface_temp      = 288;                                                   % Earth-like warmth because why not
    Planet.surface_press     = 1*ATM_2_PASCAL;                                        % Earth-like pressure because why npt      
    Planet.surface_tension   = 0.053;                                                 % liquid sulfuric acid surface tension at 288 K (DIPPR)
    Planet.kgmolwt           = CO2_atm_kgmolwt;                                       % (CO2)
    Planet.rhoa              = ideal_gas_law(Planet);                                 % Atm density from ideal gas [kg/m3] 
    Planet.nua               = vapor_CO2_dyn_viscosity_288K/Planet.rhoa;              % CO2 vapor (kinematic) viscosity [m2/s]

elseif strcmp(planet_name,'LHS-1140-b') % Cold Super-Earth
    % WATER WORLD EXOPLANET
    Planet.rho_liquid        = liquid_H2O_density_288K;                               % liquid water at 288 K (DIPPR)   
    Planet.nu_liquid         = liquid_H2O_dyn_viscosity_288K/Planet.rho_liquid;       % liquid water viscocity at 288 K (DIPPR)
    Planet.gravity           = 18.4;                                                  % Cadieux+2024a,Cadeiux+2024b (mass of 5.6 earth and 1.73 earth radius)
    Planet.surface_temp      = 288;                                                   % Earth-like warmth for liquid water on surface
    Planet.surface_press     = 0.9869*ATM_2_PASCAL;                                   % Cadeiux+2024b (1 bar for reference pressure)             
    Planet.surface_tension   = liquid_H2O_surface_tension_288K;                       % liquid water surface tension at 288 K (DIPPR)
    Planet.kgmolwt           = N2_atm_kgmolwt;                                        % (N2) (Cadeiux+2024b)
    Planet.rhoa              = ideal_gas_law(Planet);                                 % Atm density from ideal gas [kg/m3]
    Planet.nua               = vapor_N2_dyn_viscosity_288K/Planet.rhoa;               % N2 vapor (kinematic) viscosity [m2/s]

elseif strcmp(planet_name,'55-Cancri-e') % Hot Super-Earth
    
    % molten lava world a little bigger than earth
    Planet.rho_liquid       = 2450;                                                   % liquid Mount Hood Andesite at 1500 C Murase & McBirney 1973 (Fig. 10)   
    Planet.nu_liquid        = 3.162/2450;                                             % liquid Mount Hood Andesite at 1500 C Murase & McBirney 1973 (Fig. 5)
    Planet.gravity          = 22.7;                                                   % super-Earth gravity Hu+2024 (for 1.95 earth radius and 8.8 earth mass)
    Planet.surface_temp     = 1500;                                                   % a little lower than Hu+2024 predicition of 2000K but needed to use it to use Murase & McBirney 1973 values
    Planet.surface_press    = 0.9869*ATM_2_PASCAL;                                    % [1 bar] Hu+2024 (0.01 - 100 bars)     
    Planet.surface_tension  = 0.425;                                                  % liquid Mount Hood Andesite at 1500 C Nurase & McBirney 1973 (Fig. 31)
    Planet.kgmolwt          = CO2_atm_kgmolwt;                                        % (CO2)
    Planet.rhoa             = ideal_gas_law(Planet); 
    Planet.nua              = vapor_CO2_dyn_viscosity_1500K/Planet.rhoa;

elseif strcmp(planet_name,'Lake-Titicaca')
    % high altitude lake on Earth
    Planet.rho_liquid       = liquid_H2O_density_288K;                                % Liquid density [kg/m3]                               
    Planet.nu_liquid        = liquid_H2O_dyn_viscosity_288K/Planet.rho_liquid;        % liquid water viscotity at 288 K (DIPPR)                   
    Planet.gravity          = 9.81;                                                   % Gravity [m/s2]
    Planet.surface_temp     = 288;                                                    % Surface Temperature [K]
    Planet.surface_press    = 63165.48;                                               % Surface Pressure at 3812 m above sea level [Pa]       
    Planet.surface_tension  = 0.074;                                                  % Liquid Surface Tension [N/m]
    Planet.kgmolwt          = N2_atm_kgmolwt;                                         % gram molecular weight [Kgm/mol] (N2)
    Planet.rhoa             = ideal_gas_law(Planet);
    Planet.nua              = vapor_N2_dyn_viscosity_288K/Planet.rhoa;
else
    error('%s not part of default list: %s',planet_name)
end

% THRESHOLDS (m/s) (from past tests:
% EARTH                         : 2.2
% MARS-LOW                      : 1.7
% MARS-HIGH                     : 1.2
% TITAN-ONTARIOLACUS            : 0.6
% TITAN-N2                      : 0.5
% Kepler-1649b (exo-Venus)      : 5.3
% LHS-1140b (cold super Earth)  : 2.3
% 55-Cancri-e (hot super Earth) : 36.4

% e.g., wavethreshold('Earth') = 2.2


Planet.wavethreshold = containers.Map( ...
    {'Earth', ...
    'Mars-low', ...
    'Mars-high', ...
    'Titan-OntarioLacus', ...
    'Titan-N2', ...
     'Kepler-1649-b', ...
     'LHS-1140-b', ...
     '55-Cancri-e'}, ...
    [ 2.2, ...
    1.7, ...
    1.2, ...
    0.6, ...
    0.5, ...
    5.3, ...
    2.3, ...
    36.4 ] ...
);

end

