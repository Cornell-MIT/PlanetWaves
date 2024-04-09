clc
clear
close all

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETERS:
% (1) PLANET CONDITIONS
%   (a) TITAN
Titan.rho_liquid = 540;                                                    % Hydrocarbon Liquid density [kg/m3]
Titan.nu_liquid = 3e-7;                                                    % Hydrocarbon Liquid Viscocity [m2/s]
Titan.nua = 0.0126/1e4;                                                    % Titan atmospheric gas viscocity [m2/s]
Titan.gravity = 1.352;                                                     % Titan Gravity [m/s2]
Titan.surface_temp = 92;                                                   % Titan Surface Temperature [K]
Titan.surface_press = 1.5*101300;                                          % Titan Surface Pressure [Pa]
Titan.surface_tension = 0.018;                                             % Hydrocarbon Liquid Surface Tension [N/m]
Titan.name = 'Titan';

%   (b) EARTH
Earth.rho_liquid = 997;                                                    % Water Liqid Density [kg/m3]         
Earth.nu_liquid = 1.143e-6;                                                % Water Liquid Viscocity [m2/s]
Earth.nua = 1.48e-5;                                                       % Earth atmospheric gas viscocity [m2/s]
Earth.gravity = 9.81;                                                      % Earth Gravity [m/s2]
Earth.surface_temp = 288;                                                  % Earth Surface Temperature [K]
Earth.surface_press = 1*101300;                                            % Earth Surface Pressure [Pa]
Earth.surface_tension = 0.072;                                             % Water Liquid Surface Tension [N/m]
Earth.name = 'Earth';

% (2a) MODEL GEOMETRY
Model.m = 10;                                                              % Number of Grid Cells in X-Dimension
Model.n = 10;                                                              % Number of Grid Cells in Y-Dimension
Model.o = 100;                                                              % Number of Frequency bins
Model.p = 72;                                                              % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.long = 5;                                                            % longitude grid point for sampling during plotting 
Model.lat = 5;                                                             % latitude grid point for sampling during plotting
Model.gridX = 10*1000;                                                     % Grid size in X-dimension [m]
Model.gridY = 10*1000;                                                     % Grid size in Y-dimension [m]
Model.mindelt = 0.0001;                                                    % minimum time step
Model.maxdelt = 2000.0;                                                    % maximum time step
Model.time_step = 100;                                                     % Maximum Size of time step [s] -- if set too low can lead to numerical ringing
Model.num_time_steps = 10;                                                 % Length of model run (in terms of # of time steps)
Model.tolH = NaN;                                                          % tolerance threshold for maturity 
Model.cutoff_freq = 60;                                                    % cutoff frequency bin from diagnostic to advection -- if set too low can lead to numerical ringing
Model.min_freq = 0.005;                                                     % minimum frequency to model
Model.max_freq = 35;                                                       % maximum frequency to model

% (2b) BUOY SPECIFC: 
% STATION 45012:
% https://www.ndbc.noaa.gov/station_page.php?station=45012
Model.z_data = 3.6;                                                         % elevation of wind measurement [m]
deep_bathy = 237.5.*ones(Model.m,Model.n);                                  % depth of water column beneath buoy [m]
Model.bathy_map = deep_bathy;                                               % Bathymetry of model basin [m]

% (2c) TUNING PARAMETERS 
Model.tune_A1 = 0.11;                                                       % wind sea (eq. 12 Donelan+2012)
Model.tune_mss_fac = 240;
Model.tune_Sdt_fac = 0.001;
Model.tune_Sbf_fac = 0.002;
Model.tune_cotharg = 0.2;
Model.tune_n = 2.4;

% (3) NEAR-SURFACE WIND CONDITIONS
test_speeds = [3 15];                                                         % magnitude of incoming wind [m/s]
Wind.dir = 0;                                                              % direction of incoming wind [radians]

% (4) Unidirectional currents
Uniflow.East = 0;                                                          % eastward unidirectional current [m/s]
Uniflow.North = 0;                                                         % northward unidirectional current [m/s]

% (5) HOUSEKEEPING
Etc.showplots = 0;
Etc.savedata = 0;
Etc.showlog = 0;

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

planet_to_run = Earth;



for i = 1:numel(test_speeds)
	Wind.speed  = test_speeds(i);
	[myHsig(i,:),~,~,~] = makeWaves(planet_to_run,Model,Wind,Uniflow,Etc);
    if i == 1
        figure;
    end
    plot(myHsig(i,:),'-','LineWidth',3,'DisplayName', num2str(Wind.speed))
    hold on
    drawnow
end
legend('show', 'Location', 'northwest');
title(planet_to_run.name);
xlabel('time')
ylabel('wave height')

disp('run finished')

