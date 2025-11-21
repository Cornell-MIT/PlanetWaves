clc
clear
close all

% Test script for R calculations for a simple circular shoreline

addpath(genpath(fullfile('..','subroutines')));
addpath('C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\PlanetWaves\data\Titan\TitanLakes\shoreline')

%shoreline in CCW
[x,y] = make_synthetic_shoreline(20);
% remove complexity
[x,y] = make_circle(x,y,20);

% Earth
Earth.rho_s = 2.65*1000;
Earth.rho = 1.0*1000;
Earth.nu = 1e-6;
Earth.g = 9.81;
% Earth river akin to Saraswati
Earth_river.slope = 0.002;
Earth_river.width = 700;

% simple climate
wind_mag = [5 5 5];
wind_angle = [0 90 290]; % wind going from west to east = 0, positive CCW

% run model for estimated river input
[susload_dominated,bedload_dominated] = calc_riverine_flux(Earth,Earth_river);

% run model for max wave littoral transport
Qs_waves = calc_Qs_waves(x,y,wind_mag(:),wind_angle(:),Earth.rho,Earth.g);

Qs_waves(Qs_waves==0) = NaN;
R_susload = susload_dominated.Qs./Qs_waves;
R_bedload = bedload_dominated.Qs./Qs_waves;

figure;
plot(x,y,'-k')
hold on
fill(x,y,[0 0.7 0.7],'FaceAlpha',0.2)
scatter(x,y,50,Qs_waves,'filled')
colorbar
title("Qs max from waves")
axis equal padded

figure;
plot(x,y,'-k')
hold on
fill(x,y,[0 0.7 0.7],'FaceAlpha',0.2)
scatter(x,y,50,log10(R_susload),'filled')
colormap(redblue)
cb = colorbar;
ylabel(cb,'log10(R)','FontSize',16, 'Rotation',270)
clim([-max(log10(R_susload)) max(log10(R_susload))])
title("R (suspended load dominated river)")
axis equal padded

figure;
plot(x,y,'-k')
hold on
fill(x,y,[0 0.7 0.7],'FaceAlpha',0.2)
scatter(x,y,50,R_bedload,'filled')
colormap(redblue)
cb = colorbar;
ylabel(cb,'log10(R)','FontSize',16, 'Rotation',270)
clim([-max(log10(R_bedload)) max(log10(R_bedload))])
title("R (bedload dominated river)")
axis equal padded