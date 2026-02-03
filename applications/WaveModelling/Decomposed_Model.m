clc
clear
close all

% SHOW DECOMPOSITION OF MODEL 

addpath(genpath(fullfile('..','..','planetwaves')))
addpath(fullfile('..','data','Titan','TitanLakes','Bathymetries','bathtub_bathy'))
load('..\..\data\Titan\TitanLakes\Bathymetries\SAR_bathy_cleaned\ol_main_basin.mat','smoothed_ol');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
test_speeds = 3;                                                           % wind speed
time_to_run = 60*10;                                                       % time to run model
wind_direction = pi;                                                       % wind direction

[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.3);

[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);  
Etc.showplots = 1;

for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);
    Model = calc_cutoff_freq(Planet,Model,Wind);
    [avgHsig, ~, ~, ~, ~,~,~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 

 
end
%print('-vector', '-dpdf', 'Decomposed_Model')
x = 1:numel(avgHsig);
y = avgHsig;
z = zeros(size(avgHsig));
col = flip(autumn(time_to_run),1); % yellow -> red, with 61 colors (for 61 lines)
figure;
for ii = 60:60:numel(avgHsig)
    plot(x(ii),y(ii),'.','Color',col(ii,:),'MarkerSize',30)
    hold on;
end
grid on;
xlim([0 610])
ylim([0 max(y)+0.05])
xlabel('model time step [$\Delta$ t]','interpreter','latex')
ylabel('significant wave height [m]','interpreter','latex')
% print('-vector', '-dpdf', 'Decomposed_Waveheight')

figure;
pcolor(Model.bathy_map)
colorbar
title('bathymetry')
