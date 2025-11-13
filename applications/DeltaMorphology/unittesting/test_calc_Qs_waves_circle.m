clc
clear
close all

addpath(genpath(fullfile('..','subroutines')));

% shoreline in CCW
[x,y] = make_synthetic_shoreline(20);
% remove complexity, optional
[x,y] = make_circle(x,y,20);

% Earth
g = 9.81;
rho = 1000;

% simple climate
wind_mag = [5 10 5];
wind_angle = [0 90 310]; % wind going from west to east = 0, positive CCW

% run model
[a,b] = calc_Qs_waves(x,y,wind_mag(:),wind_angle(:),rho,g);


function [x,y] = make_synthetic_shoreline(N)

    make_plot = 0;

    rng(42); % random seed
    angles = linspace(0, 2*pi, N);


    radii = 5 + 0.2.*randn(1, N) + 5.*sin(3*angles);
    radii = radii.*1000;
    x = radii .* cos(angles);
    y = radii .* sin(angles);
    x = [x x(1)];
    y = [y y(1)];

    if make_plot
        figure('Name','Raw Shoreline');
        plot(x,y)
        fill(x,y,'b')
    end


end