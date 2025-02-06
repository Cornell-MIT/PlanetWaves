clc
clear
close all

% PLOT THE WAVES IN JEZERO CRATER LAKE 

addpath(fullfile('..','..','data','Mars'))
addpath(fullfile('..','..','planetwaves'))
addpath(fullfile('..','..','planetwaves','pre_analysis'))

fn = 'M20_JezeroCrater_CTXDEM_20m.tif';
lake_level = 10;

planet_to_run = 'Mars-high';
time_to_run = 60*10;
test_speeds = [1:0.5:10];
wind_direction = pi;

A = read(Tiff(fn,'r'));

% clean up Tiff
A(A<-1e10) = NaN;
A(A>-2400) = 0;
A(A>0) = 0;

% isolate crater lake
A(:,1:500) = [];
A(:,2900:end) = [];
A(1:1500,:) = [];
A(2500:end,:) = [];

grid_resolution = [20 20]; % 20 m/pix

% sz = size(A);
% X = 20.*(1:sz(2));
% Y = 20.*(1:sz(1));

buoy_loc = [580 836]; % location of Jezero delta


[A,buoy_loc,grid_resolution] = degrade_depth_resolution(A,buoy_loc,grid_resolution,0.01);

% setting water level within the crater
A(A>0) = 0;
min_original = min(min(A)); 
max_original = max(max(A));
min_new = -lake_level;
max_new = 0;
A = (A - min_original) / (max_original - min_original) * (max_new - min_new) + min_new;
A = -A;
A(A<1) = 0;


[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,A,buoy_loc);

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

save('MarsJezero.mat','myHsig','htgrid','wn','energy','cg')
make_plots(Planet,Model,Wind,test_speeds,myHsig, htgrid,energy,wn)