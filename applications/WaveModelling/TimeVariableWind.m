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
zDep = imrotate(zDep,180);

zDep_orig = zDep;
% MODEL INPUTS
planet_to_run = 'Titan-OntarioLacus';
buoy_loc = [40 35];                                                        % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
time_to_run = 60*20;                                                          % time to run model
wind_direction = 0;                                                       % wind direction

[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.3);

[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);  

make_input_map(Planet,Model,Wind)

test_speeds = make_wind_timeseries(time_to_run, [0.5, 0.5], [4 1]);

Wind.series = test_speeds;
Wind.speed = test_speeds(1);
Model = calc_cutoff_freq(Planet,Model,Wind);

[myHsig, htgrid, wn_e_spectrum, ~ , ~ , ~, ~] = makeWaveswtime(Planet, Model, Wind, Uniflow, Etc);  
if ~isempty(wn_e_spectrum{end})
    energy = squeeze(sum(wn_e_spectrum{end}.E(Model.long,Model.lat,:,:),4));
    wn = squeeze(sum(wn_e_spectrum{end}.k(Model.long,Model.lat,:,:),4));
    cg = squeeze(sum(wn_e_spectrum{end}.cg(Model.long,Model.lat,:,:),4));
end


figure
plot(myHsig)
xline(90)
yline(0.158992000000000)
yline(2.52743000000000)

%save('OntarioLacus.mat','myHsig','htgrid','wn','energy','cg')
%make_plots(Planet,Model,Wind,test_speeds,myHsig, htgrid,energy,wn)


function result = make_wind_timeseries(num_timesteps, fractions, u)

if abs(sum(fractions) - 1) > 1e-6
        error('Fractions must sum to 1.');
    end
    if length(fractions) ~= length(u)
        error('Number of fractions must match number of values.');
    end

    counts = round(fractions * num_timesteps);
    
    diff = num_timesteps - sum(counts);
    counts(end) = counts(end) + diff;

    result = [];
    for i = 1:length(counts)
        result = [result, repmat(u(i), 1, counts(i))];
    end
end
