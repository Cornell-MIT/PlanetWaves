clc
clear
close all

% Test script for R calculations for a simple circular shoreline

addpath(genpath(fullfile('..','subroutines')));
addpath('C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\PlanetWaves\data\Titan\TAMwTopo')

%shoreline in CCW
[x,y] = make_circle(1e4,10);

% Earth
Earth.rho_s = 2.65*1000;
Earth.rho = 1.0*1000;
Earth.nu = 1e-6;
Earth.g = 9.81;

Titan.rho_s = 950;
Titan.rho = 540;
Titan.nu = 6e-7;
Titan.g = 1.352;
% Earth river akin to Saraswati
Earth_river.slope = 0.002;
Earth_river.width = 700;

% simple climate
% wind_mag = 5.*ones(1,100);%[5 5 5];
% wind_angle = 90.*ones(1,100);%[0 90 290]; % wind going from west to east = 0, positive CCW

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

% run model for estimated river input
[susload_dominated,bedload_dominated] = calc_riverine_flux(Titan,Earth_river);

% run model for max wave littoral transport
Qs_waves = calc_Qs_waves(x,y,10.*wind_mag(:),wind_angle_deg(:),Titan.rho,Titan.rho_s,Titan.g);

close all

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
scatter(x,y,50,log10(R_bedload),'filled')
colormap(redblue)
cb = colorbar;
ylabel(cb,'log10(R)','FontSize',16, 'Rotation',270)
clim([-max(log10(R_bedload)) max(log10(R_bedload))])
title("R (bedload dominated river)")
axis equal padded