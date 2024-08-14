function [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_name,time_to_run,wind_dir,wind_speed,depth_profile,buoy_loc)
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

Model.tolH = NaN;                                                          % tolerance threshold for maturity
Model.min_freq = 0.05;                                                     % minimum frequency to model
Model.max_freq = 35;                                                       % maximum frequency to model
Model.z_data = 10;                                                         % elevation of wind measurement [m]
Model.bathy_map = depth_profile;
Model.num_time_steps = time_to_run;  

Wind.dir = wind_dir;
Wind.speed = wind_speed;

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
Model.mindelt = 0.0001;                                                    % seconds
Model.maxdelt = 2000.0;                                                    % seconds
  
Model.time_step = 60;                                                      % time step size [seconds]

if strcmp(planet_name,'Titan')
    % TITAN CONDITIONS
    Planet.rho_liquid = 540;                                               % Liquid density [kg/m3]
    Planet.nu_liquid = 3e-7;                                               % Liquid Viscosity [m2/s]
    Planet.nua = 0.0126/1e4;                                               % Atmospheric gas viscosity [m2/s]
    Planet.gravity = 1.352;                                                % Gravity [m/s2]
    Planet.surface_temp = 92;                                              % Surface Temperature [K]
    Planet.surface_press = 1.5*ATM_2_PASCAL;                               % Surface Pressure [Pa]
    Planet.surface_tension = 0.018;                                        % Liquid Surface Tension [N/m]
    Planet.kgmolwt = 0.028;                                                % gram molecular weight [Kgm/mol] (N2)
                                                     
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
    Planet.kgmolwt = 0.028;                                                % (N2)

elseif strcmp(planet_name,'Mars')
    % MARS CONDITIONS AT JEZERO
    Planet.rho_liquid = 997;                                                     
    Planet.nu_liquid = 1e-6;                                           
    Planet.nua = 1.4e-5;                                                  
    Planet.gravity = 3.71;                                                 
    Planet.surface_temp = 288;                                             
    Planet.surface_press = 50000;                                       
    Planet.surface_tension = 0.072;            
    Planet.kgmolwt = 0.044;                                                % (CO2)

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
    Planet.kgmolwt = 0.028;                                                % (N2)

    Model.cutoff_freq = round((15/35)*Model.Fdim);

elseif strcmp(planet_name,'55Cancrie') 
    error('not working yet')
    % molten lava world a little bigger than earth
    Planet.rho_liquid = 2450;                                                     
    Planet.nu_liquid = 12.5;                                           
    Planet.nua = 1.48e-5;                                                  
    Planet.gravity = 22.7;                                                 
    Planet.surface_temp = 1500;                                      
    Planet.surface_press = 10*100*ATM_2_PASCAL;                                       
    Planet.surface_tension = 0.44;            
    Planet.kgmolwt = 0.028;                                                % (N2)

    Model.cutoff_freq = round((15/35)*Model.Fdim);

elseif strcmp(planet_name,'C3H8')
    % propane world
    error('not working yet')
    Planet.rho_liquid = 725.30;                                                     
    Planet.nu_liquid = 5.993e-3;                                           
    Planet.nua = 2.8876e-6;                                                  
    Planet.gravity = 1.352;                                                 
    Planet.surface_temp = 92;                                      
    Planet.surface_press = 1.5*ATM_2_PASCAL;                                       
    Planet.surface_tension = 0.035923;            
    Planet.kgmolwt = 0.028;                                                % (N2)

    Model.cutoff_freq = round((12/35)*Model.Fdim);

elseif strcmp(planet_name,'N2')
    % liquid nitrogen
    Planet.rho_liquid = 734.87;                                              
    Planet.nu_liquid = 9.3244e-5;                                               
    Planet.nua = 6.432e-6;                                              
    Planet.gravity = 1.352;                                               
    Planet.surface_temp = 75;                                              
    Planet.surface_press = 1.5*ATM_2_PASCAL;                             
    Planet.surface_tension = 5.6964e-3;                                       
    Planet.kgmolwt = 0.028;                                                % (N2)                                          

    Model.cutoff_freq = round((8/35)*Model.Fdim);

elseif strcmp(planet_name,'Titan-OntarioLacus')
    % TITAN CONDITIONS at Ontario Lacus (ethane-rich)
    % 47% CH4, 40% C2H6, and 13% N2 [Mastrogiuseppe+2018] (values taken from Steckloff TitanPool for Methane Alkaline Fraction 0.47 @ 92K)
    Planet.rho_liquid = 588.15;                                               
    Planet.nu_liquid = 8.084e-7;                                               
    Planet.nua = 0.0126/1e4;                                               
    Planet.gravity = 1.352;                                               
    Planet.surface_temp = 92;                                             
    Planet.surface_press = 1.5*ATM_2_PASCAL;                             
    Planet.surface_tension = 0.032766;                                    
    Planet.kgmolwt = 0.028;                                                % (N2)

    Model.cutoff_freq = round((15/35)*Model.Fdim);

elseif strcmp(planet_name,'Titan-LigeiaMare')
    % TITAN CONDITIONS at Ligea Mare (methane-rich)
    % 69% CH4, 14% C2H6, and 17% N2 [Mastrogiuseppe+2016] (values taken from Steckloff TitanPool for Methane Alkaline Fraction 0.69 @ 92K)
    Planet.rho_liquid = 554.07;                                               
    Planet.nu_liquid = 5.8764e-7;                                              
    Planet.nua = 0.0126/1e4;                                               
    Planet.gravity = 1.352;                                               
    Planet.surface_temp = 92;                                             
    Planet.surface_press = 1.5*ATM_2_PASCAL;                             
    Planet.surface_tension = 0.028916;
    Planet.kgmolwt = 0.028;                                                % (N2)
                                  
    Model.cutoff_freq = round((15/35)*Model.Fdim);

elseif strcmp(planet_name,'Titan-CH4N2')
    % TITAN CONDITIONS: Only Methane and Nitrogen
    % [Mastrogiuseppe+2016] (values taken from Steckloff TitanPool for Methane Alkaline Fraction 1.0 @ 92K)
    Planet.rho_liquid = 510.65;                                               
    Planet.nu_liquid = 3.827e-7;                                              
    Planet.nua = 0.0126/1e4;                                               
    Planet.gravity = 1.352;                                               
    Planet.surface_temp = 92;                                             
    Planet.surface_press = 1.5*ATM_2_PASCAL;                             
    Planet.surface_tension = 0.028916;  
    Planet.kgmolwt = 0.028;                                                % (N2)
                                  
    Model.cutoff_freq = round((15/35)*Model.Fdim);
    
elseif strcmp(planet_name,'Titan-CH3H8N2')
    % TITAN CONDITIONS: Only Ethane and Nitrogen
    % [Mastrogiuseppe+2016] (values taken from Steckloff TitanPool for Methane Alkaline Fraction 0.0 @ 92K)
    Planet.rho_liquid = 654.69;                                               
    Planet.nu_liquid = 1.6501e-6;                                              
    Planet.nua = 0.0126/1e4;                                               
    Planet.gravity = 1.352;                                               
    Planet.surface_temp = 92;                                             
    Planet.surface_press = 1.5*ATM_2_PASCAL;                             
    Planet.surface_tension = 0.028916;  
    Planet.kgmolwt = 0.028;                                                % (N2)
                                  
    Model.cutoff_freq = round((15/35)*Model.Fdim);
else
    error('%s not part of default list: %s',planet_name)
end

% cutoff frequency bin from diagnostic to advection -- if set too low can lead to numerical ringing
[freqs,~] = make_frequency_vector(Model);                                  % vector of frequencies being used in wave model
for i = 1:length(Wind.speed)
    cutoff_value = (0.53*Planet.gravity)/Wind.speed(i);                    % calculate the cutoff freqency value [Hz]
    if Wind.speed(i) == 0                                                  % Check for division by zero
        Model.cutoff_freq(i) = 2;                                          % Assign 2 for undefined cutoff frequency
    else
        [~, index] = min(abs(freqs - cutoff_value));                       % Find the index of the closest frequency
        Model.cutoff_freq(i) = index;                                      % Store the index of the closest frequency
    end
end

end
