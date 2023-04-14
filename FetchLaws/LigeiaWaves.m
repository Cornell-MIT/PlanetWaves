clc
clear
close all

% make plot of signifigant wave height for Ligeia bathymetry

windspeeds = [0.4 1 3.3];

% Lorenz 2012: 78:17:5 Ethane:Methane:N2
% values for rho, nu from Titanpool at 90K

rho_lm = 630; % kg/m3
nu_lm = 1.3e-6; % m2/s

addpath('..\bathymetries')

load("Ligeia_Smoothed.mat")

bathy_map = Depth;

sigH_lm = makeWaves(windspeeds,rho_lm,nu_lm,bathy_map);


figure
plot(windspeeds,sigH_lm(:,2),'-o')
xlabel('wind speed [m/s]')
ylabel('sig H [cm]')