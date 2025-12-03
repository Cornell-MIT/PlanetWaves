clc
clear
close all

% Calculating R on Ontario Lacus

addpath(genpath(fullfile('..','subroutines')));
addpath('C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\PlanetWaves\data\Titan\TitanLakes\shoreline')
addpath('C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\PlanetWaves\data\Titan\TAMwTopo')

load('ontariolacus_shoreline.mat','lon','lati')

x = lon(2:end);
y = lati(2:end);
Titan_Radius_m = 2575*1000;
[x, y] = deg2utm(y, x, Titan_Radius_m);

% !!!! TODO make function to switch between planet conditions

% Earth
Earth.rho_s = 2.65*1000;
Earth.rho = 1.0*1000;
Earth.nu = 1e-6;
Earth.g = 9.81;

% Earth river akin to Saraswati
Earth_river.slope = 0.002;
Earth_river.width = 700;

% TAM
load('OL_winds.mat','mag_wind','angle_wind')
wind_mag = mag_wind;
wind_angle_deg = wrapTo360(round((angle_wind))); % wind going from west to east = 0, positive CCW

step_size = 100;
end_size = numel(wind_mag);

wind_mag = wind_mag(1:step_size:end_size);
wind_angle_deg = wind_angle_deg(1:step_size:end_size);

wind_mag = round(wind_mag,1);

figure;
wind_rose(wind_angle_deg,wind_mag)
title('TAM, Ontario lacus')
% 
% % run model for estimated river input
[susload_dominated,bedload_dominated] = calc_riverine_flux(Earth,Earth_river);
% 
% run model for max wave littoral transport
Qs_waves = calc_Qs_waves(x,y,wind_mag(:),wind_angle_deg(:),Earth.rho,Earth.rho_s,Earth.g);

Qs_waves(Qs_waves==0) = NaN;
R_susload = susload_dominated.Qs./Qs_waves;
R_bedload = bedload_dominated.Qs./Qs_waves;

figure;
plot(x,y,'-ok')
hold on
fill(x,y,[0 0.7 0.7],'FaceAlpha',0.2)
scatter(x,y,50,Qs_waves,'filled')
colorbar
title("Qs max from waves")
axis equal padded

figure;
plot(x,y,'-ok')
hold on
fill(x,y,[0 0.7 0.7],'FaceAlpha',0.2)
scatter(x,y,50,log10(R_susload),'filled')
colormap(redblue)
cb = colorbar;
ylabel(cb,'log10(R)','FontSize',16, 'Rotation',270)
%clim([-max(log10(R_susload)) max(log10(R_susload))])
title("R (suspended load dominated river)")
axis equal padded

figure;
plot(x,y,'-ok')
hold on
fill(x,y,[0 0.7 0.7],'FaceAlpha',0.2)
scatter(x,y,50,log10(R_bedload),'filled')
colormap(redblue)
cb = colorbar;
ylabel(cb,'log10(R)','FontSize',16, 'Rotation',270)
%clim([-max(log10(R_bedload)) max(log10(R_bedload))])
title("R (bedload dominated river)")
axis equal padded