clc
clear
close all

% unit test

% INPUT PARAMETERS:
% (1) PLANET CONDITIONS
%   (a) TITAN
Titan.rho_liquid = 465;                                                    % Hydrocarbon Liquid density [kg/m3]
Titan.nu_liquid = 0.0031/1e4;                                              % Hydrocarbon Liquid Viscocity [m2/s]
Titan.gravity = 1.352;                                                     % Titan Gravity [m/s2]
Titan.surface_temp = 92;                                                   % Titan Surface Temperature [K]
Titan.surface_press = 1.5*101300;                                          % Titan Surface Pressure [Pa]
Titan.surface_tension = 0.018;                                             % Hydrocarbon Liquid Surface Tension [N/m]


% (2) MODEL GEOMETRY
Model.m = 31;                                                              % Number of Grid Cells in X-Dimension
Model.n = 15;                                                              % Number of Grid Cells in Y-Dimension
Model.gridX = 1000.0;                                                      % Grid size in X-dimension [m]
Model.gridY = 1000.0;                                                      % Grid size in Y-dimension [m]
Model.time_step = 1;                                                       % Maximum Size of time step [s]
Model.num_time_steps = 100;                                                % Length of model run (in terms of # of time steps)
Model.bathy_map = 100.*ones(Model.m,Model.n);;                             % Bathymetry of model basin [m]

% (3) NEAR-SURFACE WIND CONDITIONS
Wind.dir = 0;                                                              % direction of incoming wind [radians]
Wind.speed = 0.4:1:3.3;                                                    % magnitude of incoming wind [m/s]

% (4) HOUSEKEEPING
Etc.showplots = 0;
Etc.savedata = 1;
Etc.showlog = 1;

% RUN MODEL
[sigH,htgrid,E_each] = makeWaves(Titan,Model,Wind,Etc); % [m]

assert()