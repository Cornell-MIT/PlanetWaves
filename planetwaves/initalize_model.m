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

ATM_2_PASCAL = 101325;

Planet.name = planet_name;

size_lake = size(depth_profile);                                    

Model.LonDim = size_lake(2);                                               % Number of Grid Cells in X-Dimension (col count)
Model.LatDim = size_lake(1);                                               % Number of Grid Cells in Y-Dimension (row count)
Model.Fdim = 50;                                                           % Number of Frequency bins
Model.Dirdim = 72;                                                         % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.long = buoy_loc(1);                                                  % longitude grid point for sampling during plotting
Model.lat = buoy_loc(2); 

Model.z_data = 10;                                                         % elevation of wind measurement [m]
Model.bathy_map = depth_profile;                                           % bathymetric map [m]

Model.tolH = NaN;                                                          % tolerance threshold for maturity
Model.min_freq = 0.05;                                                     % minimum frequency to model
Model.max_freq = 35;                                                       % maximum frequency to model


%Model.cutoff_freq = round((20/35)*Model.Fdim);                            % cutoff frequency bin from diagnostic to advection -- if set too low can lead to numerical ringing

Wind.dir = wind_dir;

Uniflow.East = 0;                                                          % eastward unidirectional current [m/s]
Uniflow.North = 0;                                                         % northward unidirectional current [m/s]

Etc.showplots = 0;
Etc.savedata = 0;

% Empirical-ish values for wave equations
Model.tune_A1 = 0.11;                                                      % wind sea (eq. 12 Donelan+2012)
Model.tune_mss_fac = 360;
Model.tune_Sdt_fac = 0.001;
Model.tune_Sbf_fac = 0.002;
Model.tune_cotharg = 0.2;
Model.tune_n = 2.4;
Model.explim = 0.1;

% sub-time step parameters for time evolution (can be optimized for particular planet condition)
Model.mindelt = 0.0001;                                                    % seconds
Model.maxdelt = 2000.0;                                                    % seconds
  
Model.time_step = 60;                                                      % time step size [seconds]
Model.num_time_steps = time_to_run;  


% VALUES FOR LIQUID DENSITY, LIQUID VISCOSITY, VAPOR (ATM) VISCOSITY, AND
% SURFACE TEMPERATURE TAKEN FROM dippr.aiche.org DATABASE
if strcmp(planet_name,'Earth')
    % EARTH CONDITIONS
    Planet.rho_liquid = 998.21;                                            % Liquid density [kg/m3]                               
    Planet.nu_liquid = 1.2e-6;                                             % Liquid Viscosity [m2/s]                        
    Planet.nua = 1.8153e-8;                                                % Atmospheric gas viscosity [m2/s] (at zero elevation: 1.81206E-5/998.21)
    Planet.gravity = 9.81;                                                 % Gravity [m/s2]
    Planet.surface_temp = 288;                                             % Surface Temperature [K]
    Planet.surface_press = 1*ATM_2_PASCAL;                                 % Surface Pressure [Pa]      
    Planet.surface_tension = 0.074;                                        % Liquid Surface Tension [N/m]
    Planet.kgmolwt = 0.028;                                                % gram molecular weight [Kgm/mol] (N2)

elseif strcmp(planet_name,'Mars-low')
    % MARS CONDITIONS AT JEZERO
    Planet.rho_liquid = 998.21;                                            % liquid water density at 288 K (DIPPR)        
    Planet.nu_liquid = 1.2e-6;                                             % liquid water viscotity at 288 K (DIPPR)
    Planet.nua = 1.4e-8;                                                   % CO2 vapor viscocity (DIPPR)         
    Planet.gravity = 3.71;                                                 % Mars gravity
    Planet.surface_temp = 288;                                             % Earth-like warmth for liquid water to exist
    Planet.surface_press = 50000;                                          % low pressure case
    Planet.surface_tension = 0.074;                                        % liquid water surface tension at 288 K (DIPPR)
    Planet.kgmolwt = 0.044;                                                % (CO2)

elseif strcmp(planet_name,'Mars-high')
    % MARS CONDITIONS AT JEZERO
    Planet.rho_liquid = 998.21;                                            % liquid water density at 288 K (DIPPR)         
    Planet.nu_liquid = 1.2e-6;                                             % liquid water viscocity at 288 K (DIPPR)
    Planet.nua = 1.4e-8;                                                   % CO2 vapor viscocity (DIPPR)
    Planet.gravity = 3.71;                                                 % Mars gravity
    Planet.surface_temp = 288;                                             % Earth-like warmth for liquid water to exist
    Planet.surface_press = 4*50000;                                        % high pressure case
    Planet.surface_tension = 0.074;                                        % liquid water surface tension at 288 K (DIPPR)
    Planet.kgmolwt = 0.044;                                                % (CO2)

elseif strcmp(planet_name,'Titan-N2')
    % TITAN CONDITIONS IN PAST/FUTURE WITH NO METHANE TO KEEP MOON WARM
    % liquid nitrogen (based on methane depletion case from Charnay+2014) for surface albedo of 0.3 with non-radiative clouds
    Planet.rho_liquid = 807.22 ;                                           % liquid N2 density at 77.5 K (DIPPR)   
    Planet.nu_liquid = 1.6e-7;                                             % liquid N2 viscocity at 77.5 (DIPPR)  
    Planet.nua = 5.4e-9;                                                   % vapor N2 viscocity at 77.5 K (DIPPR)
    Planet.gravity = 1.352;                                                % TItan gravity
    Planet.surface_temp = 77.5;                                            % Charney+2014                           
    Planet.surface_press = 0.89*ATM_2_PASCAL;                              % Charney+2014
    Planet.surface_tension = 0.0088249;                                    % liquid N2 surface tension at 77.5 K (DIPPR)   
    Planet.kgmolwt = 0.028;                                                % (N2)                                          

elseif strcmp(planet_name,'Titan-OntarioLacus')
    % TITAN CONDITIONS at Ontario Lacus (ethane-rich)
    % 47% CH4, 40% C2H6, and 13% N2 [Mastrogiuseppe+2018] (values taken from Steckloff TitanPool for Methane Alkaline Fraction 0.47 @ 92K)
    Planet.rho_liquid = 588.15;                                            % TITANPOOL   
    Planet.nu_liquid = 8.084e-7;                                           % TITANPOOL   
    Planet.nua = 6.4e-9;                                                   % vapor N2 viscocity at 92 K (DIPPR)
    Planet.gravity = 1.352;                                                % Titan gravity
    Planet.surface_temp = 92;                                              % ~avg Titan surface temp
    Planet.surface_press = 1.5*ATM_2_PASCAL;                               % ~avg Titan surface pressure
    Planet.surface_tension = 0.032766;                                     % TITANPOOL
    Planet.kgmolwt = 0.028;                                                % (N2)

elseif strcmp(planet_name,'Titan-LigeiaMare')
    % TITAN CONDITIONS at Ligea Mare (methane-rich)
    % 71% CH4, 12% C2H6, and 17% N2 [Mastrogiuseppe+2016] (values taken from Steckloff TitanPool for Methane Alkaline Fraction 0.71 @ 92K)
    Planet.rho_liquid = 551.06;                                            % TITANPOOL  
    Planet.nu_liquid = 5.685e-7;                                           % TITANPOOL   
    Planet.nua = 6.4e-9;                                                   % vapor N2 viscocity at 92 K (DIPPR)
    Planet.gravity = 1.352;                                                % Titan gravity
    Planet.surface_temp = 92;                                              % ~avg Titan surface temp
    Planet.surface_press = 1.5*ATM_2_PASCAL;                               % ~avg Titan surface pressure 
    Planet.surface_tension = 0.028363;                                     % TITANPOOL
    Planet.kgmolwt = 0.028;                                                % (N2)
                                  
elseif strcmp(planet_name,'Titan-CH4N2')
    % TITAN CONDITIONS: Only Methane and Nitrogen
    % [Mastrogiuseppe+2016] (values taken from Steckloff TitanPool for Methane Alkaline Fraction 1.0 @ 92K)
    Planet.rho_liquid = 510.65;                                            % TITANPOOL   
    Planet.nu_liquid = 3.1827e-7;                                          % TITANPOOL   
    Planet.nua = 6.4e-9;                                                   % vapor N2 viscocity at 92 K (DIPPR)
    Planet.gravity = 1.352;                                                % Titan gravity
    Planet.surface_temp = 92;                                              % ~avg Titan surface temp
    Planet.surface_press = 1.5*ATM_2_PASCAL;                               % ~avg Titan surface pressure
    Planet.surface_tension = 0.016606;                                     % TITANPOOL
    Planet.kgmolwt = 0.028;                                                % (N2)
                                      
elseif strcmp(planet_name,'Titan-CH3H8N2')
    % TITAN CONDITIONS: Only Ethane and Nitrogen
    % [Mastrogiuseppe+2016] (values taken from Steckloff TitanPool for Methane Alkaline Fraction 0.0 @ 92K)
    Planet.rho_liquid = 654.69;                                            % TITANPOOL   
    Planet.nu_liquid = 1.6501e-6;                                          % TITANPOOL    
    Planet.nua = 6.4e-9;                                                   % vapor N2 viscocity at 92 K (DIPPR)
    Planet.gravity = 1.352;                                                % Titan gravity
    Planet.surface_temp = 92;                                              % ~avg Titan surface temp
    Planet.surface_press = 1.5*ATM_2_PASCAL;                               % ~avg Titan surface pressure
    Planet.surface_tension = 0.03317;                                      % TITANPOOL
    Planet.kgmolwt = 0.028;                                                % (N2)

elseif strcmp(planet_name,'Exo-Venus') % Kepler-1649b
    % Sulfuric Acid
    Planet.rho_liquid = 1838.1;                                            % liquid sulfuric acid density at 288 K (DIPPR)                                                                         
    Planet.nu_liquid = 3.4e-5;                                             % liquid sulfuric acid viscocity at 288 K (DIPPR)
    Planet.nua = 1.4e-8;                                                   % CO2 vapor viscocity (DIPPR)
    Planet.gravity = 8.87;                                                 % Venus-like gravity
    Planet.surface_temp = 288;                                             % Earth-like warmth because why not
    Planet.surface_press = 1*ATM_2_PASCAL;                                 % Earth-like pressure because why npt      
    Planet.surface_tension = 0.053;                                        % liquid sulfuric acid surface tension at 288 K (DIPPR)
    Planet.kgmolwt = 0.044;                                                % (CO2)

elseif strcmp(planet_name,'LHS-1140b')
    % WATER WORLD EXOPLANET
    Planet.rho_liquid = 998.21;                                            % liquid water at 288 K (DIPPR)   
    Planet.nu_liquid = 1.2e-6;                                             % liquid water viscocity at 288 K (DIPPR)
    Planet.nua = 6.4e-9;                                                   % vapor N2 viscocity at 92 K (DIPPR)
    Planet.gravity = 18.4;                                                 % Cadieux+2024a,Cadeiux+2024b (mass of 5.6 earth and 1.73 earth radius)
    Planet.surface_temp = 288;                                             % Earth-like warmth for liquid water on surface
    Planet.surface_press = 0.9869*ATM_2_PASCAL;                            % Cadeiux+2024b (1 bar for reference pressure)             
    Planet.surface_tension = 0.074;                                        % liquid water surface tension at 288 K (DIPPR)
    Planet.kgmolwt = 0.028;                                                % (N2) (Cadeiux+2024b)

elseif strcmp(planet_name,'55-Cancrie') 
    % molten lava world a little bigger than earth
    Planet.rho_liquid = 2450;                                              % liquid Mount Hood Andesite at 1500 C Murase & McBirney 1973 (Fig. 10)   
    Planet.nu_liquid = 3.162/2450;                                         % liquid Mount Hood Andesite at 1500 C Murase & McBirney 1973 (Fig. 5)
    Planet.nua = 5.2e-8;                                                   % vapor CO2 at 1500 K (DIPPR)
    Planet.gravity = 22.7;                                                 % super-Earth gravity Hu+2024 (for 1.95 earth radius and 8.8 earth mass)
    Planet.surface_temp = 1500;                                            % a little lower than Hu+2024 predicition of 2000K but needed to use it to use Murase & McBirney 1973 values
    Planet.surface_press = 0.9869*ATM_2_PASCAL;                            % [1 bar] Hu+2024 (0.01 - 100 bars)     
    Planet.surface_tension = 0.425;                                        % liquid Mount Hood Andesite at 1500 C Nurase & McBirney 1973 (Fig. 31)
    Planet.kgmolwt = 0.044;                                                % (CO2)
                                  
else
    error('%s not part of default list: %s',planet_name)
end

end

