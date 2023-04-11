clc
clear
close all

% make plots of signifigant wave height for different compositions


windspeeds = 0.4:1:3.3;
bathy_map = 100.*ones(m,n);

% from titanpool at 90K

rho_methane = 540;
nu_methane = NaN;

rho_ethane = 660;
nu_ethane = NaN;

% Ontario composition of 51:38:11 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
rho_ontario = NaN;
nu_ontario = NaN;

% Punga composition of 80:0:20 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
rho_punga = NaN;
nu_punga = NaN;

sigH_methane = makeWaves(windspeeds,rho_methane,nu_methane,bathy_map);
sigH_ethane = makeWaves(windspeeds,rho_ethane,nu_ethane,bathy_map);
sigH_ontario = makeWaves(windspeeds,rho_ontario,nu_ontario,bathy_map);
sigH_punga = makeWaves(windspeeds,rho_punga,nu_punga,bathy_map);

figure
plot(u,sigH_methane,'-o')
hold on
plot(u,sigH_ethane,'-o')
plot(u,sigH_ontario,'-o')
plot(u,sigH_punga,'-o')
legend('Methane-N2','Ethane-N2','Methane-Ethane-N2 (Ontario)','Methane-Ethane-N2 (Punga)')
xlabel('wind speed [m/s]')
ylabel('sig H [cm]')