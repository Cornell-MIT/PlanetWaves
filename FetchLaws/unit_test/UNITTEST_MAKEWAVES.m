clc
clear
close all

% unit test
addpath('..')
% INPUT PARAMETERS:
% (1) PLANET CONDITIONS
%   (a) TITAN
Titan.rho_liquid = 465;                                                    % Hydrocarbon Liquid density [kg/m3]
Titan.nu_liquid = 0.0031/1e4;                                              % Hydrocarbon Liquid Viscocity [m2/s]
Titan.nua = 0.0126/1e4;                                                    % Titan atmospheric gas viscocity [m2/s]
Titan.gravity = 1.35;                                                      % Titan Gravity [m/s2]
Titan.surface_temp = 92;                                                   % Titan Surface Temperature [K]
Titan.surface_press = 1.5*101300;                                          % Titan Surface Pressure [Pa]
Titan.surface_tension = 0.018;                                             % Hydrocarbon Liquid Surface Tension [N/m]


% (2) MODEL GEOMETRY
Model.m = 31;                                                              % Number of Grid Cells in X-Dimension
Model.n = 15;                                                              % Number of Grid Cells in Y-Dimension
Model.long = 29;                                                           % longitude grid point for sampling during plotting 
Model.lat = 8;                                                             % latitude grid point for sampling during plotting
Model.o = 25;                                                              % Number of Frequency bins
Model.p = 288;                                                             % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.gridX = 1000.0;                                                      % Grid size in X-dimension [m]
Model.gridY = 1000.0;                                                      % Grid size in Y-dimension [m]
Model.mindelt = 0.0001;                                                    % minimum time step
Model.maxdelt = 2000.0;                                                    % maximum time step
Model.tolH = NaN;                                                          % tolerance threshold for maturity 
Model.time_step = 10000;                                                   % Maximum Size of time step [s]
Model.num_time_steps = 2;                                                  % Length of model run (in terms of # of time steps)
Model.bathy_map = 100.*ones(Model.m,Model.n);                              % Bathymetry of model basin [m]

% (3) NEAR-SURFACE WIND CONDITIONS
Wind.dir = 0;                                                              % direction of incoming wind [radians]
Wind.speed = 0.4:1:3.3;                                                    % magnitude of incoming wind [m/s]

% (4) Unidirectional currents
Uniflow.East = 0;                                                          % eastward unidirectional current [m/s]
Uniflow.North = 0;                                                         % northward unidirectional current [m/s]

% (4) HOUSEKEEPING
Etc.showplots = 0;
Etc.savedata = 1;
Etc.showlog = 1;

% RUN MODEL
[sigH,htgrid,E_each] = makeWaves(Titan,Model,Wind,Uniflow,Etc); % [m]

ref = runtests('refTest.m');
res = runtests('resultsTest.m');
rtref = table(ref)
rtres = table(res)