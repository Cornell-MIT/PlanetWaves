clc
clear
close all

% TEST TITAN CUTOFF FREQUENCY

addpath(fullfile('..','planetwaves'))  
addpath(fullfile('..','planetwaves','pre_analysis'))
addpath(fullfile('..','data','Titan','TitanLakes','Bathymetries','bathtub_bathy'))
load('..\data\Titan\TitanLakes\Bathymetries\bathtub_bathy\ol_bathtub_0.002000_slope.mat','zDep');

% MODEL INPUTS
planet_to_run = 'Titan';
buoy_loc = [400 800];                                                      % grid location [x,y]
grid_resolution = [10*1000 10*1000];                                       % pixel width and pixel height [m]
time_to_run = 10;                                                          % time to run model
wind_direction = 0;                                                        % wind direction
test_cutoff_freq = 10:20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL MODEL
% degrade depth profile so model doesnt take as long to run
[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.02);
% replace with uniform cube of depth equal to deepest part of Ontario
zDep = max(max(zDep)).*ones(size(zDep));
% populate model classes
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,0,zDep,buoy_loc);
% update grid resolution
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);                                               
Wind.speed = 1;                                                            % wind speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL

figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1])


for j = 1:numel(test_cutoff_freq)
    t = (j-1)/(11-1);
    mycolor = (1-t)*[1,1,0] + t*[0.5,0,0];

    Model.cutoff_freq = test_cutoff_freq(j);

    [myHsig{j}, ~,~, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  

    plot(myHsig{j},'LineWidth',2,'DisplayName',strcat('cutoff freq index =',num2str(test_cutoff_freq(j))),'Color',mycolor)
    hold on
    drawnow;
end

legend('show','Location','eastoutside')
xlabel('model time step [$\Delta$ t]','interpreter','latex')
ylabel('significant wave height [m]','interpreter','latex')
title(strcat('Cut-off Frequency on ',Planet.name))
