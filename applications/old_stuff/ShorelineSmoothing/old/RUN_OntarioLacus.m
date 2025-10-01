clc
clear
close all
% 
% % SHORELINE SMOOTHING ON ONTARIO LACUS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % HOUSEKEEPING
addpath(fullfile('..','..','..','planetaryrivermeander\external\scripts\wavelet_Perron'))
addpath(fullfile('..','..','..','planetaryrivermeander\external\scripts\wavelet_Torrence'))
addpath(fullfile('..','..','..','planetaryrivermeander\src\user_function'))

addpath(fullfile('..','..','planetwaves'))  
addpath(fullfile('..','..','planetwaves','pre_analysis'))
addpath(fullfile('..','..','planetwaves','post_analysis'))
addpath(fullfile('..','..','data','Titan','TitanLakes','Bathymetries','bathtub_bathy'))

% simple bathtub model of depths
load('ontario_lacus_bathtub_0.002.mat','x','y','Xmesh','Ymesh','Depth')

gap = 700;
[x,y] = smooth_path(x,y,2);
[x,y] = even_spacing(x,y,gap);

% A = [x y];
% AA = sqrt( sum( abs( diff( A ) ).^2, 2 ) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP BATHYMETRY FOR WAVE MODEL


% MODEL INPUTS
planet_to_run = 'Titan-OntarioLacus';
buoy_loc = [750 150];                                                      % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
test_speeds = 3;                                                           % wind speed
time_to_run = 60;                                                          % time to run model
wind_direction = pi;                                                        % wind direction

figure;
subplot(1, 2, 1);
surf(Xmesh, Ymesh, Depth,'EdgeColor','none');
view(2)
hold on;
plot3(x, y, max(Depth(:)) * ones(size(x)), 'r', 'LineWidth', 2);
title('Original Depth and Path');
xlabel('X'); ylabel('Y'); zlabel('Depth');

[zDep, buoy_loc, new_resolution, Xmesh, Ymesh, ~, ~] = degrade_depth_mesh(Depth, buoy_loc, grid_resolution, 1/50, Xmesh, Ymesh, x, y);

% Display results
% zDep(zDep<0) = NaN;

buoy_loc = [7 15];
zDep = round(zDep);
zDep(zDep==0) = NaN;


[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction(1),zDep,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);  

subplot(1, 2, 2);
surf(Xmesh, Ymesh, zDep);
hold on;
plot3(x, y, max(zDep(:)) * ones(size(x)), 'r', 'LineWidth', 2);
title('Resized Depth and Path');
xlabel('X'); ylabel('Y'); zlabel('Depth');
view(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRIDDED PLOT FOR DEPTH

WIDTH = (max(x)-min(x))/Model.LonDim;                                      % Width of each grid cell (# cells * width of 1 cell = width of all cells)
HEIGHT = (max(y)-min(y))/Model.LatDim;                                     % Height of each grid cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN WAVE MODEL (need H, T) OF MAIN ENERGY COMPONENT
Wind.speed = test_speeds;
Wind.dir = wind_direction;

Model = calc_cutoff_freq(Planet,Model,Wind);

make_input_map(Planet,Model,Wind)

[~,~,~,~,~,~,PeakWave] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  

% EXTRACT GRID OF WAVE HEIGHTS AND PERIODS
sig_wave = invert_attribute(PeakWave);

POI = 1;
h = plot_wave_grid(x,y,WIDTH,HEIGHT,Model,sig_wave);
hold on
plot(x(POI),y(POI),'.r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE Qs and Psi for each sub-grid point from wave grid
[wind_dir,wave_front_angle,shoreline_angle,relative_angle] = calc_shoreline_angle(x,y,Wind);

wave_energy = calc_wave_energy(x,y,sig_wave,WIDTH,HEIGHT,Model,Planet);
[Qs,mydiff] = shoreline_smoothing(x,y,sig_wave,relative_angle,WIDTH,HEIGHT,Model);
sinu = calc_sinuosity(x,y,200);


figure;
scatter(x,y,100,sinu,"filled")
colorbar
title('sinuosity')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE POINTS NORMAL RELATIVE TO THE DIRECTION OF THE WIND

my_angle = point_rel2_wind(x,y,wind_dir);

% jet_wrap = vertcat(jet,flipud(jet));
% figure;
% 
% scatter(x,y,100,rad2deg(my_angle),"filled")
% colormap(jet_wrap);
% colorbar

for i = 1:numel(x)
    if my_angle(i) > pi/2 && my_angle(i) < 3*pi/2
        along_wind(i) = false;
    else
        along_wind(i) = true;
    end
end


alongwind_diff = mydiff(along_wind);
alongwind_sin = sinu(along_wind);
not_alongwind_diff = mydiff(~along_wind);
not_alongwind_sin = sinu(~along_wind);
alongwind_wave_energy = wave_energy(along_wind);
not_alongwind_wave_energy = wave_energy(~along_wind);

figure;
scatter(x,y,100,mydiff,'ok')
hold on
scatter(x(along_wind),y(along_wind),100,alongwind_diff,"filled")
colorbar

[~,dans_along,xans_along,yans_along] = kde2d([alongwind_diff alongwind_sin']);

[~,dans_notalong,xans_notalong,yans_notalong] = kde2d([not_alongwind_diff not_alongwind_sin']);
figure;
contour3(xans_notalong,yans_notalong,dans_notalong,50,'-r')
view(2)
hold on;
contour3(xans_along,yans_along,dans_along,50,'-k')

xlabel('diffusivity')
ylabel('sinuoisty')
ylim([0 1])



figure;
scatter(x,y,100,wave_energy,'ok')
hold on
scatter(x(along_wind),y(along_wind),100,alongwind_wave_energy./max(alongwind_wave_energy),"filled")
colorbar


[~,dans_along,xans_along,yans_along] = kde2d([alongwind_wave_energy alongwind_sin']);

[~,dans_notalong,xans_notalong,yans_notalong] = kde2d([not_alongwind_wave_energy not_alongwind_sin']);
figure;
contour3(xans_notalong,yans_notalong,dans_notalong,50,'-r')
view(2)
hold on;
contour3(xans_along,yans_along,dans_along,50,'-k')

xlabel('wave energy')
ylabel('sinuoisty')
ylim([0 1])


function newval = closest_power2(val)
    newval = pow2(floor(log2(val)));
end

function [xEven, yEven] = even_spacing(x, y,gap)

    
    total_length = arclength(x,y,'linear');

    p = interparc(0:(gap/total_length):1,x,y,'linear');
    xEven = p(:,1);
    yEven = p(:,2);
    
end




