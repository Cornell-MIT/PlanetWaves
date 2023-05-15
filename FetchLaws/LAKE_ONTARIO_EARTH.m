clc
clear
close all

% Lake Ontario Bathymetry: https://www.ncei.noaa.gov/products/great-lakes-bathymetry
% Citation: National Geophysical Data Center, 1999. Bathymetry of Lake Ontario. National Geophyiscal Data Center, NOAA. https://doi.org/10.7289/V56H4FBH
addpath('..\bathymetries\')
% export DEM values into array for UMWM bathymetry
D = imread('OntarioDEM1.tif');
D(D==D(1,1)) = 0; % set land values to 0 (depth is positive)
% figure;
% contourf(D,20)
% xlabel('x')
% ylabel('y')
% colormap(hot)
% xlim([0 249])
% ylim([0 70])

[m,n] = size(D);

windspeeds = 7.7; % 15 knots (waves build to 0.6-1.2 meters)

rho_water = 997;
nu_water = 1e-6;

wind_dir = 0;


[sigH,E_spec] = makeWaves(windspeeds,wind_dir,rho_water,nu_water,D,1,100); % [m]
