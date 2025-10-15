clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Housekeeping

% Ontario Lacus shoreline
addpath(fullfile('..','..','\data\Titan\TitanLakes\shoreline'))
% Juan's winds
addpath(fullfile('..','..','data/Titan/TAMwTopo/')) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate max flux from waves

% Load in shoreline
load('ontariolacus_shoreline.mat')
x = lon; y = lati;
x(1) = []; y(1) = [];
%[x, y] = make_circle(x,y); % make a simple shape for testing

% Calculate the regional shoreline orientaiton at each point
window_size_ang = 5; % size of window to define the regional shoreline orientation over
shoreline_angle = calc_regional_shoreline_angle(x, y,window_size_ang);


% Specify wave climate
wave_mean = 180; % degrees CCW from E
wave_std = 500;
wave = 0:359;
E_pdf = vonMises(deg2rad(wave_mean),deg2rad(wave_std),deg2rad(wave));

% Calculate wave flux of sediment
wave_sed_flux = calc_wave_flux(E_pdf,rad2deg(shoreline_angle)); % this is not specific to a size?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate max flux from largest unobservable rivers

Titan_cold.rho_s = 0.95*1000; 
Titan_cold.rho = 0.67*1000;
Titan_cold.nu = 0.003/10000;
Titan_cold.g = 1.35;

vidflumina.coldbedload_D50 = [3.8 10]/100;
vidflumina.susload_D50 = 6.35e-5;
vidflumina.width = [100 175];
vidflumina.slope = [0.0011 0.0015];

[river_susload,river_bedload] = riverine_flux(Titan_cold.rho_s,Titan_cold.rho,Titan_cold.nu,Titan_cold.g,vidflumina.coldbedload_D50(1),vidflumina.susload_D50,vidflumina.width(1),NaN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate fluvial dominance R (Qriver/Qwaves)

R = river_bedload.Qs./wave_sed_flux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTS 
% Plot of shoreline with quiver showing the window size calculating the
% regional angle over
figure('Name','window size for angle estimation');
plot(x,y,'-k','MarkerFaceColor','k')
hold on
quiver(x(1),y(1),x(window_size_ang) - x(1),y(window_size_ang) - y(1),'k')
title('window size')
axis equal padded

% Plot showing the regional shoreline angle at each point
figure('Name','shoreline angle')
scatter(x,y,10,(rad2deg(shoreline_angle)))
colorbar;
title('regional shoreline angle')
axis equal padded

% Plot showing the maximum flux from the waves
figure('Name','Qsmax')
scatter(x,y,10,wave_sed_flux)
colormap(jet)
colorbar
title('Qsmax of waves')
axis equal padded


viridis_top = viridis;
viridis_top = viridis_top(1:end/2,:);
plasma_bottom = plasma;
plasma_bottom = plasma_bottom(end/2:end,:);
distinct_cmap = [viridis_top;plasma_bottom];

% Plot showing fluvial dominance ratio R 
% Large R -> River Dominated
% Small R -> Wave Dominated
% R = 1 -> Waves and rivers equal
figure('Name','Delta Morph R')
scatter(x,y,10,log10(R))
colormap(distinct_cmap)
colorbar;
clim([-10 10])
title('log10 R')
axis equal padded





