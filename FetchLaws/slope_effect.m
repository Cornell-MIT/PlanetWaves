clc
clear
close all

% make plots of signifigant wave height for slopes

windspeeds = [1 2 3];

m = 31;                                                                    % number of gridpoints in x-direction
n = 15;                                                                    % number of gridpoints in y-direction

bathy_map = 100.*ones(m,n);

% from titanpool at 90K

rho_methane = 540;
nu_methane = 3e-7; % m2/s



sigH_flat = makeWaves(windspeeds,rho_methane,nu_methane,bathy_map,1,1000); % [m]


figure;
plot(1:length(sigH_flat),sigH_flat(1,:),'-o')
hold on;
plot(1:length(sigH_flat),sigH_flat(2,:),'-o')
plot(1:length(sigH_flat),sigH_flat(3,:),'-o')
legend('u = 1 m/s','u = 2 m/s','u = 3 m/s')




