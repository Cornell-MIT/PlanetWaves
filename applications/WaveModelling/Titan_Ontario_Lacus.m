clc
clear
close all

% PLOT WAVES IN ONTARIO LACUS

addpath(fullfile('..','..','planetwaves'))  
addpath(fullfile('..','..','planetwaves','pre_analysis'))
addpath(fullfile('..','..','data','Titan','TitanLakes','Bathymetries','bathtub_bathy'))
load('..\..\data\Titan\TitanLakes\Bathymetries\SAR_bathy_cleaned\ol_main_basin.mat','smoothed_ol');

% isolate main basin of interest
zDep = smoothed_ol;
zDep = imrotate(zDep,180);
zDep = imrotate(zDep,-90);
zDep(:,80:end) = [];
zDep(1:95,:) = [];
zDep(90:end,:) = [];
zDep = imrotate(zDep,-90);

zDep_orig = zDep;
% MODEL INPUTS
planet_to_run = 'Titan-OntarioLacus';
buoy_loc = [60 55];                                                        % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
test_speeds = [0.5:0.5:4];                                                % wind speed
time_to_run = 60*10;                                                          % time to run model
wind_direction = pi;                                                       % wind direction

[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.3);

[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);  

make_input_map(Planet,Model,Wind)


for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);
    Model = calc_cutoff_freq(Planet,Model,Wind);

    [myHsig{i}, htgrid{i}, wn_e_spectrum, ~ , ~ , ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
    if ~isempty(wn_e_spectrum{end})
        energy{i} = squeeze(sum(wn_e_spectrum{end}.E(Model.long,Model.lat,:,:),4));
        wn{i} = squeeze(sum(wn_e_spectrum{end}.k(Model.long,Model.lat,:,:),4));
        cg{i} = squeeze(sum(wn_e_spectrum{end}.cg(Model.long,Model.lat,:,:),4));
    end
end
save('OntarioLacus.mat','myHsig','htgrid','wn','energy','cg')
make_plots(Planet,Model,Wind,test_speeds,myHsig, htgrid,energy,wn)