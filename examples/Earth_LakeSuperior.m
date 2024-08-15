clc
clear
close all

% WAVES IN LAKE SUPERIOR CONDITIONS (BUT UNIFORM DEPTH)

warning('need to streamline process of moving outputs of find_fetch.py to matlab scripts to run model')

addpath(fullfile('..','planetwaves'))  
addpath(fullfile('..','planetwaves','pre_analysis'))
addpath(fullfile('..','data','Earth','GreatLakes','LakeSuperior'))
% from find_fetch.py 
load('..\data\Earth\GreatLakes\LakeSuperior\BathyData\LakeSuperior_cleaned.mat')
zDep = -squeeze(LS);

% MODEL INPUTS
planet_to_run = 'Earth';
buoy_loc = [6618 1729];                                                    % grid location [x,y]
% from find_fetch.py
grid_resolution = [4542.948547909539 92.66280063299297];                   % pixel width and pixel height [m]
test_speeds = 0:20;                                                        % wind speed
time_to_run = 60*10;                                                       % time to run model
wind_direction = 0;                                                        % wind direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL MODEL
% degrade depth profile so model doesnt take as long to run
[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.02);
zDep = max(max(zDep)).*ones(12,12);
% populate model classes
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
% adjust anenometer height to buoy
Model.z_data = 3.6;

% update grid resolution
Model.gridX = 20*1000;                                              
Model.gridY = 20*1000;                                               

Model.long = 6;                                                  
Model.lat = 6; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL

% Preallocate cell arrays to store results
myHsig = cell(1, numel(test_speeds));
htgrid = cell(1, numel(test_speeds));
E_spec = cell(1, numel(test_speeds));
Cg = cell(1,numel(test_speeds));

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PLOT RESULTS
save('lakesuperior_run.mat','myHsig','energy','wn','cg')

make_plots(Planet,Model,Wind,test_speeds,myHsig, htgrid,energy,wn)