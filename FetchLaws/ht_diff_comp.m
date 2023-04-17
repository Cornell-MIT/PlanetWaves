clc
clear
close all

% make plots of signifigant wave height for different compositions


<<<<<<< Updated upstream
windspeeds = [0.4 1];
=======
windspeeds = [0.4:.4:3.4];
>>>>>>> Stashed changes

m = 31;                                                                    % number of gridpoints in x-direction
n = 15;                                                                    % number of gridpoints in y-direction

bathy_map = 100.*ones(m,n);

<<<<<<< Updated upstream
% from titanpool at 90K
=======
% % from titanpool at 90K
>>>>>>> Stashed changes

rho_methane = 540;
nu_methane = 3e-7; % m2/s

rho_ethane = 660; % Titanpool and Hayes 2012
nu_ethane = (0.0011)/rho_ethane; % kinematic viscocity (titanpool Hayes 2012)

<<<<<<< Updated upstream
=======
%from titanpool at 92K (steckloff et al., 2020)
comp = [0:10:100];
rho_min = 654.091;
rho_max = 512.692;
rho_comp = (((rho_max-rho_min)./(max(comp)-min(comp))).*comp) + rho_min;
nu_min = nu_ethane;
nu_max = nu_methane;
nu_comp = (((nu_max-nu_min)./(max(comp)-min(comp))).*comp) + nu_min;

>>>>>>> Stashed changes
% % Ontario composition of 51:38:11 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
% rho_ontario = NaN;
% nu_ontario = NaN;
% 
% % Punga composition of 80:0:20 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
% rho_punga = NaN;
% nu_punga = NaN;

<<<<<<< Updated upstream
sigH_methane = makeWaves(windspeeds,rho_methane,nu_methane,bathy_map);
sigH_ethane = makeWaves(windspeeds,rho_ethane,nu_ethane,bathy_map);
% sigH_ontario = makeWaves(windspeeds,rho_ontario,nu_ontario,bathy_map);
% sigH_punga = makeWaves(windspeeds,rho_punga,nu_punga,bathy_map);

figure
plot(windspeeds,sigH_methane(:,2),'-o')
hold on
plot(windspeeds,sigH_ethane(:,2),'-o')
legend('Methane-N2','Ethane-N2')
xlabel('wind speed [m/s]')
ylabel('sig H [cm]')
=======
% sigH_methane = makeWaves(windspeeds,rho_methane,nu_methane,bathy_map);
% sigH_ethane = makeWaves(windspeeds,rho_ethane,nu_ethane,bathy_map);
% sigH_ontario = makeWaves(windspeeds,rho_ontario,nu_ontario,bathy_map);
% sigH_punga = makeWaves(windspeeds,rho_punga,nu_punga,bathy_map);
sigH_100ethane_0methane = makeWaves(windspeeds,rho_comp(:,1),nu_comp(:,1),bathy_map);
sigH_90ethane_10methane = makeWaves(windspeeds,rho_comp(:,2),nu_comp(:,2),bathy_map);
sigH_80ethane_20methane = makeWaves(windspeeds,rho_comp(:,3),nu_comp(:,3),bathy_map);
sigH_70ethane_30methane = makeWaves(windspeeds,rho_comp(:,4),nu_comp(:,4),bathy_map);
sigH_60ethane_40methane = makeWaves(windspeeds,rho_comp(:,5),nu_comp(:,5),bathy_map);
sigH_50ethane_50methane = makeWaves(windspeeds,rho_comp(:,6),nu_comp(:,6),bathy_map);
sigH_40ethane_60methane = makeWaves(windspeeds,rho_comp(:,7),nu_comp(:,7),bathy_map);
sigH_30ethane_70methane = makeWaves(windspeeds,rho_comp(:,8),nu_comp(:,8),bathy_map);
sigH_20ethane_80methane = makeWaves(windspeeds,rho_comp(:,9),nu_comp(:,9),bathy_map);
sigH_10ethane_90methane = makeWaves(windspeeds,rho_comp(:,10),nu_comp(:,10),bathy_map);
sigH_0ethane_100methane = makeWaves(windspeeds,rho_comp(:,11),nu_comp(:,11),bathy_map);

% figure
% plot(windspeeds,sigH_methane(:,2),'-o')
% hold on
% plot(windspeeds,sigH_ethane(:,2),'-o')
% legend('Methane-N2','Ethane-N2')
% xlabel('wind speed [m/s]')
% ylabel('sig H [cm]')
% hold off

figure
plot(windspeeds,sigH_100ethane_0methane(:,2),'-o', 'Color', [66/255 0/255 83/255] )
hold on
plot(windspeeds,sigH_90ethane_10methane(:,2),'-o', 'Color', [69/255 32/255 118/255])
plot(windspeeds,sigH_80ethane_20methane(:,2),'-o', 'Color', [62/255 66/255 137/255])
plot(windspeeds,sigH_70ethane_30methane(:,2),'-o', 'Color', [49/255 96/255 143/255])
plot(windspeeds,sigH_60ethane_40methane(:,2),'-o', 'Color', [38/255 120/255 144/255])
plot(windspeeds,sigH_50ethane_50methane(:,2),'-o', 'Color', [26/255 145/255 141/255])
plot(windspeeds,sigH_40ethane_60methane(:,2),'-o', 'Color', [30/255 169/255 133/255])
plot(windspeeds,sigH_30ethane_70methane(:,2),'-o', 'Color', [66/255 192/255 114/255])
plot(windspeeds,sigH_20ethane_80methane(:,2),'-o', 'Color', [123/255 210/255 79/255])
plot(windspeeds,sigH_10ethane_90methane(:,2),'-o', 'Color', [187/255 224/255 35/255])
plot(windspeeds,sigH_0ethane_100methane(:,2),'-o', 'Color', [253/255 232/255 33/255])
legend('100% Ethane - 0% Methane','90% Ethane - 10% Methane','80% Ethane - 20% Methane', ...
    '70% Ethane - 30% Methane','60% Ethane - 40% Methane','50% Ethane - 50% Methane', ...
    '40% Ethane - 60% Methane','30% Ethane - 70% Methane','20% Ethane - 80% Methane', ...
    '10% Ethane - 90% Methane','0% Ethane - 100% Methane')
xlabel('wind speed [m/s]')
ylabel('sig H [cm]')
title('windspeed vs sig H at 92K')
hold off


>>>>>>> Stashed changes
