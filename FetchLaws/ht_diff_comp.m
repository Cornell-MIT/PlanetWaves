clc
clear
close all

% make plots of signifigant wave height for different compositions


windspeeds = [0.5 1];

m = 31;                                                                    % number of gridpoints in x-direction
n = 15;                                                                    % number of gridpoints in y-direction

bathy_map = 100.*ones(m,n);

% from titanpool at 90K

rho_methane = 540;
nu_methane = 3e-7; % m2/s

rho_ethane = 660; % Titanpool and Hayes 2012
nu_ethane = (0.0011)/rho_ethane; % kinematic viscocity (titanpool Hayes 2012)

% % Ontario composition of 51:38:11 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
% rho_ontario = NaN;
% nu_ontario = NaN;
% 
% % Punga composition of 80:0:20 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
% rho_punga = NaN;
% nu_punga = NaN;

sigH_methane = makeWaves(windspeeds,rho_methane,nu_methane,bathy_map);
% sigH_ethane = makeWaves(windspeeds,rho_ethane,nu_ethane,bathy_map);
% sigH_ontario = makeWaves(windspeeds,rho_ontario,nu_ontario,bathy_map);
% sigH_punga = makeWaves(windspeeds,rho_punga,nu_punga,bathy_map);

figure
plot(windspeeds,sigH_methane(:,1),'-o')
hold on
plot(windspeeds,sigH_methane(:,2),'-o')
legend('T1','T2')
xlabel('wind speed [m/s]')
ylabel('sig H [m]')

