clc
clear
close all

addpath(genpath(fullfile('..','subroutines')));

%  square shoreline w equal sides
x = [0 10 10 0 0]*1000; % meters
y = [0 0 10 10 0]*1000;

% Earth
g = 9.81;
rho = 1000;

% simple climate 
wind_mag   = [5 5 10];    % m/s
wind_angle = [0 90 45];   % where wind blows to, deg CCW from east

% run model
[Qs_total, Qs_time] = calc_Qs_waves(x, y, wind_mag, wind_angle, rho, g);
