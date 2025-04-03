clc
clear
close all

% making waves within a basin of uniform depth everywhere instead of sloping
addpath('C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\umwm_titan\data\Titan\TitanLakes\Bathymetries\bathtub_bathy')
addpath('C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\umwm_titan\applications\ShorelineSmoothing')
% WAVES FOR 8 DIRECTIONS
wind_direction = 0:45:315;
wind_direction = deg2rad(wind_direction);
fne = {'0deg.mat','45deg.mat','90deg.mat','135deg.mat','180deg.mat','225deg.mat','270deg.mat','315deg.mat'};

load('asylake1.mat','x','y');

[Xmesh,Ymesh,zDep] = make_bathtub_lake(1e-3,[x y]);
close all

zDep(~isnan(zDep)) = 80;

planet_to_run = 'Titan-OntarioLacus';

buoy_loc = [500 500];
grid_resolution = [10000 10000];
test_speeds = 3;                                                           % wind speed
time_to_run = 60;                                                          % time to run model

[zDep, buoy_loc, new_resolution, Xmesh, Ymesh, ~, ~] = degrade_depth_mesh(zDep, buoy_loc, grid_resolution, 1/20, Xmesh, Ymesh, x, y);

zDep(zDep==0) = 0;
zDep(zDep~=0) = 80;

for i = 1:numel(fne)

    [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction(i),zDep,buoy_loc);
    
    Model.gridX = grid_resolution(1);                                              
    Model.gridY = grid_resolution(2);  
    
    WIDTH = (max(x)-min(x))/Model.LonDim;                                      % Width of each grid cell (# cells * width of 1 cell = width of all cells)
    HEIGHT = (max(y)-min(y))/Model.LatDim;                                     % Height of each grid cell
    
    Wind.speed = test_speeds;
    Wind.dir = wind_direction(i);
    Model = calc_cutoff_freq(Planet,Model,Wind);
    
    make_input_map(Planet,Model,Wind)
    
    [~,~,~,~,~,~,PeakWave] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
    
    % EXTRACT GRID OF WAVE HEIGHTS AND PERIODS
    sig_wave = invert_attribute(PeakWave);
    save(convertCharsToStrings(fne{i}),'sig_wave')
end

