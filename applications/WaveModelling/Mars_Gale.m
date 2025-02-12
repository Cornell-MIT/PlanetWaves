clc
clear
close all

addpath(fullfile('..','..','data','Mars','StevensRubin2022'))
% addpath(fullfile('..','..','planetwaves'))

%downwind
x = linspace(0, 20000, 500); 
% crosswind
y = 1:3; 

depth = zeros(length(y), length(x));

% slope is discontinous
fetch = 20*1000;
slope_break = 18800;
flat_region = x <= slope_break; 
slope_region = x > slope_break;

% depth only varying in downwind
depth(:, flat_region) = -10; % Constant depth of 10m

% Linear slope in the last 1.2 km 
depth(:, slope_region) = repmat(-10 * (1 - (x(slope_region) - slope_break) / (fetch-slope_break)), length(y), 1);

% Create meshgrid for plotting
[X, Y] = meshgrid(x / 1000, y); % Convert to km

% Plot surface
figure;
h = pcolor(X, Y, depth);
set(h,'EdgeColor','none');
xlabel('downwind (km)');
ylabel('crosswind (km)');
zlabel('Depth (m)');
title('Gale (Rubin+2022)');
colorbar;
hold on;
xline(0:1.2:18)


load('stevensrubins2022.mat')
x_SWAN = (marswaveparams.distance_km);
H_SWAN = (marswaveparams.hsig_m);

x_SWAN_p6_u7p5 = x_SWAN(marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==7.5);
x_SWAN_p6_u15 = x_SWAN(marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==15);
x_SWAN_p6_u30 = x_SWAN(marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==30);
x_SWAN_p6_u60 = x_SWAN(marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==7.5);
x_SWAN_p6_u120 = x_SWAN(marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==7.5);
H_SWAN_p6_u7p5 = H_SWAN(marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==7.5);
H_SWAN_p6_u15 = H_SWAN(marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==15);
H_SWAN_p6_u30 = H_SWAN(marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==30);
H_SWAN_p6_u60 = H_SWAN(marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==60);
H_SWAN_p6_u120 = H_SWAN(marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==120);

figure;
plot(x_SWAN_p6_u7p5,H_SWAN_p6_u7p5,'--r');
hold on
plot()

% downwind_sz = 1.2*1000;
% num_downwind_cells = numel(0:1.2:20);
% crosswind_sz = 10*1000;
% num_crosswind_cells = 10;
% lake_level = 10;
% 
% gale_crater = lake_level*ones(num_crosswind_cells,num_downwind_cells);
% figure;
% pcolor(1:num_downwind_cells,1:num_crosswind_cells,gale_crater)
% 
% planet_to_run = 'Mars-high';
% time_to_run = 60;
% test_speeds = [3];
% wind_direction = 0;
% 
% buoy_loc = [6 6];
% [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,gale_crater,buoy_loc);
% Model.gridX = downwind_sz;                                              
% Model.gridY = 10*1000;  
% Model.Fdim = 24;
% Model.min_freq = 0.05;                                                     % minimum frequency to model
% Model.max_freq = 1;  
% Planet.surface_press = 60*1000;
% Planet.surface_temp = 273;
% make_input_map(Planet,Model,Wind)
% 
% for i = 1:numel(test_speeds)
% 
%     Wind.speed = test_speeds(i);
%     Model = calc_cutoff_freq(Planet,Model,Wind);
% 
%     [~, htgrid{i},~, ~ , ~ , ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
% 
% end
% 
