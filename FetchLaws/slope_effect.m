clc
clear
close all

% make plots of signifigant wave height for slopes

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

% (2) MODEL GEOMETRY
Model.m = 10;                                                              % Number of Grid Cells in X-Dimension
Model.n = 10;                                                              % Number of Grid Cells in Y-Dimension
Model.o = 25;                                                              % Number of Frequency bins
Model.p = 288;                                                             % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.long = 5;                                                            % longitude grid point for sampling during plotting 
Model.lat = 5;                                                             % latitude grid point for sampling during plotting
Model.gridX = 1000;                                                        % Grid size in X-dimension [m]
Model.gridY = 1000;                                                        % Grid size in Y-dimension [m]
Model.mindelt = 0.0001;                                                    % minimum time step
Model.maxdelt = 2000.0;                                                    % maximum time step
Model.time_step = 10;                                                      % Maximum Size of time step [s]
Model.num_time_steps = 10;                                                 % Length of model run (in terms of # of time steps)
Model.tolH = NaN;                                                          % tolerance threshold for maturity 

% define a bathymetry with a constant slope
deep_bathy = 100.*ones(Model.m,Model.n);
Model.bathy_map = deep_bathy;                                              % Bathymetry of model basin [m]

% (3) NEAR-SURFACE WIND CONDITIONS
Wind.dir = 0;                                                              % direction of incoming wind [radians]
Wind.speed = 0:0.1:2;                                                    % magnitude of incoming wind [m/s]

% (4) Unidirectional currents
Uniflow.East = 0;                                                          % eastward unidirectional current [m/s]
Uniflow.North = 0;                                                         % northward unidirectional current [m/s]

% (5) HOUSEKEEPING
Etc.showplots = 0;
Etc.savedata = 0;
Etc.showlog = 1;

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

planet_to_run = Titan;

% RUN MODEL
[DsigH,Dhtgrid,DE_each] = makeWaves(planet_to_run,Model,Wind,Uniflow,Etc); 

PM_H = 0.22.*(Wind.speed.^2)./(planet_to_run.gravity);                             % Pierson-Moskowitz Sig Wave Height

figure;
plot(Wind.speed,DsigH(:,end),'-k','LineWidth',5)
hold on;
plot(Wind.speed,PM_H,'--r','LineWidth',5)
grid on;
xlabel('Wind Speed [m/s]')
ylabel('sig H [m]')
legend('Titan-UMWM','Pierson-Moskowitz')
title(planet_to_run.name)

figure;
[~,c] = contourf(1:Model.num_time_steps,Wind.speed,DsigH,'showtext',1);
colormap(autumn)
c.LineWidth = 3;
cb = colorbar;
cb.Label.String = 'sig H [m]';
xlabel('# Time Steps (1000 s)')
ylabel('Wind Speed [m/s]')
title(planet_to_run.name)


