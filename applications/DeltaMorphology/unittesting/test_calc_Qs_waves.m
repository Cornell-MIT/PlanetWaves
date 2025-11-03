clc
clear
close all

addpath(genpath(fullfile('..','subroutines')));



[x,y] = make_synthetic_shoreline(20);
%[x,y] = make_complicated_shoreline(100);
g = 9.81;
rho = 1000;

wind_mag = [5 10];
wind_angle = 0:10:359; % winding going from west to east = 0, positive CCW

[a,b] = calc_Qs_waves(x,y,wind_mag(:),wind_angle(:),rho,g);


function [x,y] = make_synthetic_shoreline(N)

    rng(42); % random seed
    angles = linspace(0, 2*pi, N);

 
    radii = 5 + 0.2.*randn(1, N) + 5.*sin(3*angles);
    radii = radii.*1000;
    x = radii .* cos(angles);
    y = radii .* sin(angles);
    x = [x x(1)];
    y = [y y(1)];

    figure('Name','Raw Shoreline');
    plot(x,y)
    fill(x,y,'b')

end