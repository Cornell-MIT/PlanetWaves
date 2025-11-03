
clc
clear
close all


% Time steps
T = 100;
% Synthetic waves 
%wind_mag = 5*rand(1,T);
wind_mag = 5.*ones(1,T);
wind_angle = 180*rand(1,T);

figure;
tiledlayout('horizontal')
nexttile
plot(1:numel(wind_angle),wind_angle)
xlabel('wind direction')
nexttile
plot(1:numel(wind_mag),wind_mag)
xlabel('wind mag')

% synthetic shoreline (closed polygon)
N = 200;
x = [10 30 50 60 55 40 25 15 10];  % x-coordinates
y = [10 15 10 25 40 45 35 30 10];  % y-coordinates
dist = [0, cumsum(sqrt(diff(x).^2 + diff(y).^2))];
xi = linspace(0, dist(end), N);
yi = spline(dist, y, xi);
xi = spline(dist, x, xi);
x = xi;
y = yi;

x = circshift(x,-100);
y = circshift(y,-100);

% calculate summed wave flux
[Qsmax_total, Qsmax_ts] = calc_wave_flux_per_wave(wind_angle, wind_mag, x, y);

disp('Qsmax_total:');
disp(Qsmax_total);

disp('Qsmax_ts:');
disp(Qsmax_ts);


figure;
plot(1:numel(x), Qsmax_total, '-o');
xlabel('Shoreline point index');
ylabel('Qsmax total');
title('Total Qsmax per shoreline');

figure;
imagesc(1:T, 1:numel(x), Qsmax_ts);
xlabel('Time step');
ylabel('Shoreline point index');
colorbar;
title('Time series of Qsmax per shoreline');
