clc
clear
close all

% checking model against Banfield+2015 wave tank experiments
addpath(fullfile('..','planetwaves'))  
addpath(fullfile('..','data','Mars'))

% TABLE 1
WaveTank = readtable('Banfield2015_table1.xlsx'); % [atm_pressure wind_speed sigH]

% TABLE 2
OceanA.rho_liquid = 1010;
OceanA.nu_liquid = 0.0018/OceanA.rho_liquid;
OceanA.surface_tension = 0.0722;
OceanB.rho_liquid = 1180;
OceanB.nu_liquid = 0.0220/OceanA.rho_liquid;
OceanB.surface_tension = 0.0810;
OceanC.rho_liquid = 1360;
OceanC.nu_liquid = 0.0420/OceanA.rho_liquid;
OceanC.surface_tension = 0.0900;

figure('units','normalized','outerposition',[0 0 1 1])
scatter(WaveTank.wind_speed_m_s,WaveTank.sig_H_mm./1000,50,WaveTank.atm_pressure_mbar.*1000,'filled')
grid on;
title('Banfield+2015 Wave Tank Measurements for different atmospheres [bar] (table 1)')
xlabel('wind speed m/s')
ylabel('sigH m')
colormap(winter)
colorbar


% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETERS:
% (1) PLANET CONDITIONS
%   (b) EARTH
Earth.rho_liquid = 997;                                                    % Water Liqid Density [kg/m3]       
Earth.nu_liquid = 1.143e-6;                                                % Water Liquid Viscosity [m2/s]
Earth.nua = 1.48e-5;                                                       % Earth atmospheric gas viscosity [m2/s]
Earth.gravity = 9.81;                                                      % Earth Gravity [m/s2]
Earth.surface_temp = 288;                                                  % Earth Surface Temperature [K]
Earth.surface_press = 1*101300;                                            % Earth Surface Pressure [Pa]
Earth.surface_tension = 0.072;                                             % Water Liquid Surface Tension [N/m]
Earth.name = 'Earth';
%   (c) Paleo-Mars
Mars = Earth;                      
Mars.gravity = 3.71;
Mars.surface_temp = 288;
Mars.surface_pressure = 1*101300;
Mars.name = 'Mars';
% (2a) MODEL GEOMETRY
Model.LonDim = 10;                                                         % Number of Grid Cells in X-Dimension (col count)
Model.LatDim = 10;                                                         % Number of Grid Cells in Y-Dimension (row count)
Model.Fdim = 35;                                                           % Number of Frequency bins
Model.Dirdim = 72;                                                         % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.long = 6;                                                            % longitude grid point for sampling during plotting
Model.lat = 6;                                                             % latitude grid point for sampling during plotting
Model.gridX = 1;                                                           % Grid size in X-dimension [m]
Model.gridY = 1;                                                           % Grid size in Y-dimension [m]
Model.mindelt = 0.0001;                                                    % minimum time step
Model.maxdelt = 1;                                                    % maximum time step
Model.time_step = 0.1;                                                     % Maximum Size of time step [s] -- if set too low can lead to numerical ringing
Model.num_time_steps = 200;                                                % Length of model run (in terms of # of time steps)
Model.tolH = NaN;                                                          % tolerance threshold for maturity
Model.cutoff_freq = 15;                                                    % cutoff frequency bin from diagnostic to advection -- if set too low can lead to numerical ringing
Model.min_freq = 0.05;                                                     % minimum frequency to model
Model.max_freq = 35;                                                       % maximum frequency to model
% (2b) BUOY SPECIFC:
% STATION 45012:
% https://www.ndbc.noaa.gov/station_page.php?station=45012
Model.z_data = (28.5/100);                                                  % elevation of wind measurement [m]
Model.bathy_map = (16/100).*ones(Model.LonDim,Model.LatDim);
% (2c) TUNING PARAMETERS
Model.tune_A1 = 0.11;                                                      % wind sea (eq. 12 Donelan+2012)
Model.tune_mss_fac = 360;
Model.tune_Sdt_fac = 0.001;
Model.tune_Sbf_fac = 0.002;
Model.tune_cotharg = 0.2;
Model.tune_n = 2.4;
% (3) NEAR-SURFACE WIND CONDITIONS
test_speeds = 3;%sort(unique(WaveTank.wind_speed_m_s'));                             % magnitude of incoming wind [m/s]
Wind.dir = 0;                                                              % direction of incoming wind [radians]
% (4) Unidirectional currents
Uniflow.East = 0;                                                          % eastward unidirectional current [m/s]
Uniflow.North = 0;                                                         % northward unidirectional current [m/s]
% (5) HOUSEKEEPING
Etc.showplots = 0;
Etc.savedata = 0;
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
figure; 
imagesc(Model.bathy_map)
colormap cool; 
hold on
contour(Model.bathy_map,'--k')
plot(Model.long,Model.lat,'or','MarkerFaceColor','r')
colorbar

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% RUN THE MODEL
planet_to_run = Earth;

% Preallocate cell arrays to store results
myHsig = cell(1, numel(test_speeds));
htgrid = cell(1, numel(test_speeds));
E_spec = cell(1, numel(test_speeds));
myHsig = cell(1, numel(test_speeds));
htgrid = cell(1, numel(test_speeds));
E_spec = cell(1, numel(test_speeds));

warning('need to loop over atm pressure, wind speed, and temperature. fill this in and parallelize')

parfor i = 1:numel(test_speeds)
    % Create a local copy of the Wind variable for each iteration so can be parallelized on local machine
    Wind_local = Wind;
    Wind_local.speed = test_speeds(i);
    
    [myHsig{i}, htgrid{i}, E_spec{i}, ~] = makeWaves(planet_to_run, Model, Wind_local, Uniflow, Etc);   % run model
    
end
disp('all winds completed ')

% plot results

figure
for i = 1:numel(test_speeds)
    plot(myHsig{i}, '-', 'LineWidth', 3, 'DisplayName', num2str(test_speeds(i)))
    hold on
end
grid on;
legend('show', 'Location', 'northwest','interpreter','latex');
title(['Waves on',' ',planet_to_run.name],'interpreter','latex');
xlabel('model time step [$\Delta$ t]','interpreter','latex')
ylabel('significant wave height [m]','interpreter','latex')





