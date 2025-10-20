clc
clear
close all

total = 430652;
duration_years = 10;
%stopi = round(total/(10/duration_years));
stopi = total;
% ONTARIO LACUS
loni = 34;
lati = 3; 
% LIGEIA MARE
% loni = 45; 
% lati = 31; 


% 1 Titan Year = 29.5 Earth Year = 258420 hours
% 1 Titan Day = 15.95 Earth days = 382.8 hours
% 1 Measurement per 6 hours

% Measurements per Day = (382.8 / 6) = 63.8
% Measurements per Year = (258420 / 6) = 43070


fn = "C:\Users\Owner\OneDrive\Documents\00_Main\Work\MIT\External_Data\TAM.10mwinds.L24.table950.k1e-4.y181-190.6-hourly.nc";

fileinfo = ncinfo(fn);
lon = ncread(fn, 'lon');
lat = ncread(fn, 'lat');
%time = ncread(fn, 'time');
uc = ncread(fn, 'ucomp',[1 1 1],[length(lon) length(lat) stopi]);
vc = ncread(fn, 'vcomp',[1 1 1],[length(lon) length(lat) stopi]);



u_ol = squeeze(uc(loni,lati,:));
v_ol = squeeze(vc(loni,lati,:));

mag_wind = sqrt(u_ol.^2 + v_ol.^2);
angle_wind = mod(rad2deg(atan2(v_ol,u_ol)) + 360,360);


figure('Name','Wave Energy PDF');
plot((1:length(mag_wind))./64,mag_wind);
xlabel('days')
ylabel('|u|')

figure
wind_rose(angle_wind,mag_wind)


% Parameters
Fs = 64; %samples_per_day 
days_per_year = 10759;
N= length(mag_wind); 

% Compute the FFT

Y = fft(mag_wind);

% Plot the single-sided normalized amplitude spectrum
figure;
ax = gca();

plot(1./(Fs/N*(0:N-1)),abs(Y),'-r','LineWidth',1.5);

xline(1/2, ':k', 'LineWidth', 1, 'Label', '1/2 Day');
xline(0.5*days_per_year, ':k', 'LineWidth', 1, 'Label', '1/2 Year');

title(['(lon,lat): (', num2str(lon(loni)), ',' num2str(lat(lati)),')'])
xlabel('1/freq (day)');
ylabel("|fft(X)|")
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
