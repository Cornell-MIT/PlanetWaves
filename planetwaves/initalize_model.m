function [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_name,time_to_run,wind_dir,depth_profile,buoy_loc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function will populate model parameters with default values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT:
%       planet_name     :   name of planet [string] e.g. = 'Titan'
%       time_to_run     :   number of time steps [#]
%       wind_dir        :   direction of wind [radians]
%       depth_profile   :   2D array of depths [m]
%       buoy_loc        :   grid location (x,y) of buoy for measurement
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

Model.cutoff_freq = round((20/35)*Model.Fdim);                             % cutoff frequency bin from diagnostic to advection -- if set too low can lead to numerical ringing

Model.tolH = NaN;                                                          % tolerance threshold for maturity
Model.min_freq = 0.05;                                                     % minimum frequency to model
Model.max_freq = 35;                                                       % maximum frequency to model
Model.z_data = 10;                                                         % elevation of wind measurement [m]
Model.bathy_map = depth_profile;
Model.num_time_steps = time_to_run;  

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

% sub-time step parameters for time evolution (can be optimized for particular planet condition)
Model.mindelt = 0.0001;                                               
Model.maxdelt = 2000.0;                                                   
Model.time_step = 100;       

if strcmp(planet_name,'Titan')
    % TITAN CONDITIONS
    Planet.rho_liquid = 540;                                               % Hydrocarbon Liquid density [kg/m3]
    Planet.nu_liquid = 3e-7;                                               % Hydrocarbon Liquid Viscosity [m2/s]
    Planet.nua = 0.0126/1e4;                                               % Titan atmospheric gas viscosity [m2/s]
    Planet.gravity = 1.352;                                                % Titan Gravity [m/s2]
    Planet.surface_temp = 92;                                              % Titan Surface Temperature [K]
    Planet.surface_press = 1.5*ATM_2_PASCAL;                               % Titan Surface Pressure [Pa]
    Planet.surface_tension = 0.018;                                        % Hydrocarbon Liquid Surface Tension [N/m]

    Model.time_step = 50;                                                  % Maximum Size of time step [s] -- if set too low can lead to numerical ringing
   Model.cutoff_freq = round((15/35)*Model.Fdim);
   
elseif strcmp(planet_name,'Earth')
    % EARTH CONDITIONS
    Planet.rho_liquid = 997;                                                     
    Planet.nu_liquid = 1.143e-6;                                           
    Planet.nua = 1.48e-5;                                                  
    Planet.gravity = 9.81;                                                 
    Planet.surface_temp = 288;                                             
    Planet.surface_press = 1*ATM_2_PASCAL;                                       
    Planet.surface_tension = 0.072;            
    
elseif strcmp(planet_name,'Mars')
    % MARS CONDITIONS AT JEZERO
    Planet.rho_liquid = 997;                                                     
    Planet.nu_liquid = 1e-6;                                           
    Planet.nua = 1.4e-5;                                                  
    Planet.gravity = 3.71;                                                 
    Planet.surface_temp = 288;                                             
    Planet.surface_press = 50000;                                       
    Planet.surface_tension = 0.072;            
                                                                                               
    Model.time_step = 50;
    Model.cutoff_freq = round((15/35)*Model.Fdim);

elseif strcmp(planet_name,'Exo-Venus') 
    % Sulfuric Acid
    Planet.rho_liquid = 790.03;                                                     
    Planet.nu_liquid = 1.26e-4;                                           
    Planet.nua = 1.24e-5;                                                  
    Planet.gravity = 8.87;                                                 
    Planet.surface_temp = 293.15;                                             
    Planet.surface_press = 1*ATM_2_PASCAL;                                       
    Planet.surface_tension = 0.012322;            
    % 
    % Model.mindelt = 0.001;                                                 % minimum time step
    % Model.time_step = 50;
    Model.cutoff_freq = round((15/35)*Model.Fdim);
else
    error('%s not part of default list: %s',planet_name)
end




end