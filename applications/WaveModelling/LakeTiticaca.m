clc
clear
close all

% WAVES IN LAKE TITICACA (UNIFORM DEPTH)
% rho_air = 0.785 kg/m3;
% avg_depth = 107 m
% max_width =  190 km
% max_length = 80 km 

addpath(fullfile('..','..','planetwaves'))  
addpath(fullfile('..','..','planetwaves','pre_analysis'))


% MODEL INPUTS
planet_to_run = 'Lake-Titicaca';
                                                 
test_speeds = 1:10;                                                        % wind speed
time_to_run = 60*10;                                                       % time to run model
wind_direction = 0;                                                        % wind direction

zDep = 107.*ones(19,8); 
buoy_loc = [4 9];
grid_resolution = [10*1000 10*1000];         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL MODEL
% populate model classes
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);                                               

Model.long = buoy_loc(1);                                                  
Model.lat = buoy_loc(2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL

% Preallocate cell arrays to store results
myHsig = cell(1, numel(test_speeds));
htgrid = cell(1, numel(test_speeds));
E_spec = cell(1, numel(test_speeds));
Cg = cell(1,numel(test_speeds));

figure;
hold on;
for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);
    Model = calc_cutoff_freq(Planet,Model,Wind);

    [myHsig{i}, htgrid{i}, wn_e_spectrum, ~ , ~ , ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
    if ~isempty(wn_e_spectrum{end})
        energy{i} = squeeze(sum(wn_e_spectrum{end}.E(Model.long,Model.lat,:,:),4));
        wn{i} = squeeze(sum(wn_e_spectrum{end}.k(Model.long,Model.lat,:,:),4));
        cg{i} = squeeze(sum(wn_e_spectrum{end}.cg(Model.long,Model.lat,:,:),4));
    end

    wv_ht(i) = myHsig{i}(end);
    plot(test_speeds(i),wv_ht(i),'rs')
end

plot(test_speeds,wv_ht,'-rs')
planet_to_run = 'Earth';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL MODEL
% populate model classes
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);                                               

Model.long = buoy_loc(1);                                                  
Model.lat = buoy_loc(2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL

for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);
    Model = calc_cutoff_freq(Planet,Model,Wind);

    [myHsig_Earth{i}, ~, ~, ~ , ~ , ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
    
    wv_ht_normal(i) = myHsig_Earth{i}(end);
    plot(test_speeds(i),wv_ht_normal(i),'ks')
end
plot(test_speeds,wv_ht_normal,'-ks')