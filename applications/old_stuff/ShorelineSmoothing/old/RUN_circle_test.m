clc
clear
close all

% SHORELINE SMOOTHING TESTED ON A CIRCLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOUSEKEEPING
addpath(fullfile('..','..','planetwaves'))  
addpath(fullfile('..','..','planetwaves','pre_analysis'))
addpath(fullfile('..','..','planetwaves','post_analysis'))
addpath(fullfile('..','..','data','Titan','TitanLakes','Bathymetries','bathtub_bathy'))
% simple bathtub model of depths
load('..\..\data\Titan\TitanLakes\Bathymetries\bathtub_bathy\ol_bathtub_0.002000_slope.mat','zDep')
% hi-res coordinates of shoreline
load('..\..\data\Titan\TitanLakes\shoreline\OL_SHORELINE.mat','X_cor','Y_cor')
x = X_cor;y = Y_cor;
[x,y] = smooth_path(x,y,10);
[x,y] = make_circle(x,y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP BATHYMETRY FOR WAVE MODEL

zDep = imrotate(zDep,180);

% MODEL INPUTS
planet_to_run = 'Titan-OntarioLacus';
buoy_loc = [630 255];                                                      % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
test_speeds = 3;                                                         % wind speed
time_to_run = 10;                                                         % time to run model
wind_direction = pi;                                                     % wind direction

[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.02);

zDep(zDep<1) = NaN;
[a,b] = size(zDep);
zDep = max(max(zDep)).*ones(a,b);

[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction(1),zDep,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);  

make_input_map(Planet,Model,Wind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRIDDED PLOT FOR DEPTH

WIDTH = (max(x)-min(x))/Model.LonDim;                                             % Width of each grid cell (# cells * width of 1 cell = width of all cells)
HEIGHT = (max(y)-min(y))/Model.LatDim;                                             % Height of each grid cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN WAVE MODEL (need H, T) OF MAIN ENERGY COMPONENT
Wind.speed = test_speeds;
Wind.dir = wind_direction;

Model = calc_cutoff_freq(Planet,Model,Wind);

[~, ~, ~, ~ , ~ , ~, PeakWave] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  

% EXTRACT GRID OF WAVE HEIGHTS AND PERIODS
sig_wave = PeakWave';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
POI = 51;
% CALCULATE Qs and Psi for each sub-grid point from wave grid
[wind_dir,wave_front_angle,shoreline_angle,relative_angle] = calc_shoreline_angle(x,y,Wind);

[Qs,psi] = shoreline_smoothing(x,y,sig_wave,relative_angle,WIDTH,HEIGHT,Model);
sin = calc_sinuosity(x,y,5);

figure;
scatter(x,y,50,'k')
hold on;
scatter(x,y,50,psi,'filled')
plot(x(POI),y(POI),'.r')

fprintf('wind: %f\n',rad2deg(wind_dir));
fprintf('wave: %f\n',rad2deg(wave_front_angle(POI)))
fprintf('shore: %f\n',rad2deg(shoreline_angle(POI)))
fprintf('relative angle: %f\n',rad2deg(relative_angle(POI)));
fprintf('dif: %f\n',psi(POI))

figure; scatter(x,y,50,rad2deg(relative_angle),'filled'); colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PLOT STUFF
% f1 = plot_wave_grid(x,y,WIDTH,HEIGHT,Model);
% plot(x(~isnan(theta)),y(~isnan(theta)),'or','MarkerFaceColor','r')
% plot(x(isnan(theta)),y(isnan(theta)),'ok','MarkerFaceColor','k')
% 
% scatter(x,y,100,Qs,"filled")
% colorbar
% title('Qs')
% 
% f2 = plot_wave_grid(x,y,WIDTH,HEIGHT,Model);
% scatter(x,y,100,psi,"filled")
% colorbar
% title('\Psi')
% 
% 
% figure;
% scatter(x,y,100,sin,"filled")
% colorbar
% title('sinuosity')
% 
% figure;
% plot(psi,sin,'ok','MarkerFaceColor','k')
% xlabel('\Psi')
% ylabel('sin')
% 
% 
% 
