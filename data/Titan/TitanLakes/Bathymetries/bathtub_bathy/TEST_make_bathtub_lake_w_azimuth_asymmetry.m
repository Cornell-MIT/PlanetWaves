clc
clear
close all

addpath(fullfile('..','..','shoreline')) 

load('OL_SHORELINE.mat','X_cor','Y_cor')

min_slope = 0.46e-3;
max_slope = 4.85e-3;

x = X_cor;
x(1) = [];
y = Y_cor;
y(1) = [];
shoreline = [x y];


bath_slope1 = linspace(min_slope,max_slope,round(numel(x)/2));
bath_slope2 = linspace(max_slope,min_slope,round(numel(x)/2));
bath_slope = [bath_slope1 bath_slope2];

[Xmesh, Ymesh, zDep] = make_bathtub_lake_w_azimuth_asymmetry(bath_slope, shoreline);

set(gca,'XDir','reverse')
set(gca,'YDir','reverse')