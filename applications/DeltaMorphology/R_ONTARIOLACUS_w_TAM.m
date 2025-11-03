clc
clear
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Housekeeping

% Ontario Lacus shoreline
addpath(fullfile('..','..','\data\Titan\TitanLakes\shoreline'))

% Juan's GCM winds
fl = 'C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\PlanetWaves\data\Titan\TAMwTopo';
load('OL_winds.mat','mag_wind','angle_wind')

time_slice = 1:100;
mag_wind = mag_wind(time_slice);
angle_wind = angle_wind(time_slice);

angle_wind = mod(angle_wind+180,360); % angle wind is coming from
angle_wind_rad = deg2rad(angle_wind);
angle_wind_deg = angle_wind;
clearvars angle_wind

figure;
wind_rose(angle_wind_deg,mag_wind)
figure;
histogram(mag_wind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate max flux from waves

% Load in shoreline
load('ontariolacus_shoreline.mat')
x = lon; y = lati;
x(1) = []; y(1) = [];
[x, y] = make_circle(x,y); % make a simple shape for testing



% wave_flux_sum output will have dimensions M for shoreline
% wave_flux_per_time output will have dimensions M for shoreline and N for time
[wave_flux_sum, wave_flux_dt] = calc_wave_flux_per_wave(angle_wind_deg, mag_wind, x,y);

figure
scatter(x,y,50,wave_flux_sum)
colorbar