clc
clear
close all

load('ls_low_res.mat','lake_superior_low_res')
imagesc(lake_superior_low_res)

% SPLIT INTO JOB ARRAY FOR HPC, WIND SPEEDS RUN INDEPENDANT

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
Earth.surface_temp = 283;                                                  % Earth Surface Temperature [K]
Earth.surface_press = 1*102000;                                            % Earth Surface Pressure [Pa]
Earth.surface_tension = 0.072;                                             % Water Liquid Surface Tension [N/m]
Earth.name = 'Earth';

% (2a) MODEL GEOMETRY
Model.m = 50;                                                              % Number of Grid Cells in X-Dimension
Model.n = 10;                                                             % Number of Grid Cells in Y-Dimension
Model.o = 35;                                                              % Number of Frequency bins
Model.p = 40;                                                              % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.long = 5;                                                            % longitude grid point for sampling during plotting 
Model.lat = 25;                                                            % latitude grid point for sampling during plotting
Model.gridX = 1000;                                                        % Grid size in X-dimension [m]
Model.gridY = 1000;                                                        % Grid size in Y-dimension [m]
Model.mindelt = 0.0001;                                                    % minimum time step
Model.maxdelt = 2000.0;                                                    % maximum time step
Model.time_step = 1000;                                                      % Maximum Size of time step [s]
Model.num_time_steps = 10;                                                 % Length of model run (in terms of # of time steps)
Model.tolH = NaN;                                                          % tolerance threshold for maturity 
Model.cutoff_freq = 15;                                                    % cutoff frequency bin from diagnostic to advection
Model.min_freq = 0.05;                                                     % minimum frequency to model
Model.max_freq = 35;                                                       % maximum frequency to model

% (2b) BUOY SPECIFC: 
% STATION 45012:
% https://www.ndbc.noaa.gov/station_page.php?station=45012
Model.z_data = 10;                                                        % elevation of wind measurement [m]
Model.bathy_map = zDep;                                              % Bathymetry of model basin [m] 

% (2c) TUNING PARAMETERS 
Model.tune_A1 = 0.11;
Model.tune_mss_fac = 2000;
Model.tune_Sdt_fac = 0.001;
Model.tune_Sbf_fac = 0.002;
Model.tune_cotharg = 0.2;
Model.tune_n = 2.4;

% (3) NEAR-SURFACE WIND CONDITIONS
test_speeds = [1 2 3];                                                  % magnitude of incoming wind [m/s]
Wind.dir = 0;                                                              % direction of incoming wind [radians]

% (4) Unidirectional currents
Uniflow.East = 0;                                                          % eastward unidirectional current [m/s]
Uniflow.North = 0;                                                         % northward unidirectional current [m/s]

% (5) HOUSEKEEPING
Etc.showplots = 0;
Etc.savedata = 0;
Etc.showlog = 0;

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


planet_to_run = Titan;

%% PARALLELIZE PROCESS 
% >> eval(pRun('testWaves_cluster',Np,'grid')
% In serial >> eval(pRun('testWaves_cluster',1,'grid')
% With 2 processors in parallel >> eval(pRun('testWaves_cluster',1,'grid')

% SPLIT INTO JOB ARRAY 
PARALLEL = true; % flag to toggle parallelism

% make map
map1 = 1; % for serial case

if PARALLEL
% in Hsig matrix the columns are time steps and rows are wind speeds (e.g. 3 speeds for 10 steps is 3x10 array where
% the rows are independent of one another so can be distributed among Np proceessors [Np 1]
       map1 = map([Np 1],{},0:Np-1); % {} = block processor, Np = # of processor when job submitted
end
% apply map
Hsig = zeros(numel(test_speeds),Model.num_time_steps,map1); % distributed matrix DMAT

%retrieve local part of global index
iglobal = global_ind(Hsig,1);
% retrieve local portion of distributed matrix
myHsig = local(Hsig); % local cmd copies local piece of distributed matric to processor

% take in taskID (PID) as the first argument 
%my_task_id = varargin{1}; % type = double
% take in number of tasks as the second argument
%num_tasks = varargin{2}; % type = double
% split up original list of wind speeds between different tasks
%my_winds = test_speeds(my_task_id:num_tasks:length(test_speeds));

%% RUN MODEL
for i_local = 1:length(iglobal)
       % determine global index of local iteration	
	iglob = iglobal(i_local);
	Wind.speed  = test_speeds(iglob);
	[myHsig(i_local,:),~,~] = makeWaves(planet_to_run,Model,Wind,Uniflow,Etc);
end	
% gather local arrays from each processor
% store local portion in global matrix
Hsig = put_local(Hsig,myHsig);
% gather data to leader procesor
Hsig_final = agg(Hsig);

save('OntarioBathtub.mat','Hsig_final')

disp('Run finished')


