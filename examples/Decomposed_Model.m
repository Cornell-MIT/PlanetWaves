clc
clear
close all

% MODELED VS OBSERVED AT LAKE SUPERIOR DEEPWATER BUOY

addpath(fullfile('..','planetwaves'))  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL
test_speeds = 3;
planet_to_run = 'Earth';
time_to_run = 10;   
wind_direction = 0;      
grid_resolution = [20*1000 20*1000];
zDep = 273.5.*ones(12,12);
buoy_loc = [6 6];
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,0,zDep,buoy_loc);
Model.z_data = 3.6;
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);
Model.tolH = NaN;
Etc.showplots = 1;

for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);

    [avgHsig, ~, E_spec, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 

 
end


