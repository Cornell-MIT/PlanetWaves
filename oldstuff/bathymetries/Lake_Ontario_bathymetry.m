clc
clear
close all
% Lake Ontario Bathymetry: https://www.ncei.noaa.gov/products/great-lakes-bathymetry
% Citation: National Geophysical Data Center, 1999. Bathymetry of Lake Ontario. National Geophyiscal Data Center, NOAA. https://doi.org/10.7289/V56H4FBH

% export DEM values into array for UMWM bathymetry
X = imread('OntarioDEM1.tif');
X(X==X(1,1)) = 0; % set land values to 0 (depth is positive)
figure;
contourf(X,20)
xlabel('x')
ylabel('y')
colormap(hot)
xlim([0 249])
ylim([0 70])