clc
clear
close all

% make plots of signifigant wave height for different compositions


windspeeds = [0.4 1 3.3 10];


m = 31;                                                                    % number of gridpoints in x-direction
n = 15;                                                                    % number of gridpoints in y-direction

bathy_map = 100.*ones(m,n);

%% from titanpool at 90K (220 min for 2 windspeeds; 440 min for 4 windspeeds)

rho_methane = 540;
nu_methane = 3e-7; % m2/s

rho_ethane = 660; % Titanpool and Hayes 2012
nu_ethane = (0.0011)/rho_ethane; % kinematic viscocity (titanpool Hayes 2012)

% %% methane:ethane:nitrogen from titanpool at 92K (steckloff et al., 2020) and
% %NIST
% 
% % 0:96:4.27
% rho_0_96_4 = 652.735;
% sfct_0_96_4 = 1.209E-3;
% nu_0_96_4 = 1.209E-3;
% sigH_0_96_4 = makeWaves(windspeeds,rho_0_96_4,sfct_0_96_4,nu_0_96_4,bathy_map);
% 
% % 10:86:4.96
% rho_10_86_4 = 641.081;
% sfct_10_86_4 = 3.562E-2;
% nu_10_86_4 = 1.106E-3;
% sigH_10_86_4 = makeWaves(windspeeds,rho_10_86_4,sfct_10_86_4,nu_10_86_4,bathy_map);
% 
% % 19:76:5.53
% rho_19_76_5 = 626.643;
% sfct_19_76_5 = 3.353E-2;
% nu_19_76_5 = 1.004E-3;
% sigH_19_76_5 = makeWaves(windspeeds,rho_19_76_5,sfct_19_76_5,nu_19_76_5,bathy_map);
% 
% % 28:66:6.35
% rho_28_66_6 = 612.61;
% sfct_28_66_6 = 3.144E-2;
% nu_28_66_6 = 9.02E-4;
% sigH_28_66_6 = makeWaves(windspeeds,rho_28_66_6,sfct_28_66_6,nu_28_66_6,bathy_map);
% 
% % 37:55:7.61
% rho_37_55_8 = 597.556;
% sfct_37_55_8 = 2.935E-2;
% nu_37_55_8 = 7.998E-4;
% sigH_37_55_8 = makeWaves(windspeeds,rho_37_55_8,sfct_37_55_8,nu_37_55_8,bathy_map);
% 
% % 45:45:9.42
% rho_45_45_10 = 582.749;
% sfct_45_45_10 = 2.726E-2;
% nu_45_45_10 = 6.976E-4;
% sigH_45_45_10 = makeWaves(windspeeds,rho_45_45_10,sfct_45_45_10,nu_45_45_10,bathy_map);
% 
% % 53:35:11.52
% rho_53_35_12 = 567.498;
% sfct_53_35_12 = 2.517E-2;
% nu_53_35_12 = 5.954E-4;
% sigH_53_35_12 = makeWaves(windspeeds,rho_53_35_12,sfct_53_35_12,nu_53_35_12,bathy_map);
% 
% % 60:26:14.27
% rho_60_26_14 = 552.24;
% sfct_60_26_14 = 2.308E-2;
% nu_60_26_14 = 4.932E-4;
% sigH_60_26_14 = makeWaves(windspeeds,rho_60_26_14,sfct_60_26_14,nu_60_26_14,bathy_map);
% 
% % 66:16:17.70
% rho_66_16_18 = 538.596;
% sfct_66_16_18 = 2.098E-2;
% nu_66_16_18 = 3.910E-4;
% sigH_66_16_18 = makeWaves(windspeeds,rho_66_16_18,sfct_66_16_18,nu_66_16_18,bathy_map);
% 
% % 71:8:21.59
% rho_71_8_21 = 524.596;
% sfct_71_8_21 = 1.889E-2;
% nu_71_8_21 = 1.106E-3;
% sigH_71_8_21 = makeWaves(windspeeds,rho_71_8_21,sfct_71_8_21,nu_71_8_21,bathy_map);
% 
% % 75:0:25.24
% rho_75_0_25 = 512.577;
% sfct_75_0_25 = 1.68E-2;
% nu_75_0_25 = 1.867E-4;
% sigH_75_0_25 = makeWaves(windspeeds,rho_75_0_25,sfct_75_0_25,nu_75_0_25,bathy_map);

%% lake compositions from NIST and sources (80 min for 2 windspeeds; 160 min for 4 windspeeds)
% % Ontario composition of 51:38:11 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
rho_ontario = 476.962;
sfct_ontario = 2.579E-2;
nu_ontario = 6.261E-4;
sigH_ontario = makeWaves(windspeeds,rho_ontario,sfct_ontario,nu_ontario,bathy_map);

% % Punga composition of 80:0:20 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
rho_punga = 361.006;
sfct_punga = 1.680E-2;
nu_punga = 1.866E-4;
sigH_punga = makeWaves(windspeeds,rho_punga,sfct_punga,nu_punga,bathy_map);

% % Ligia max methane composition of 100:0:0 percent
% methane:ethane:nitrogen from LeGall 2016
rho_ligiamax = 449.83;
sfct_ligiamax = 1.68E-2;
nu_ligiamax = 1.866E-4;
sigH_ligiamax = makeWaves(windspeeds,rho_ligiamax,sfct_ligiamax,nu_ligiamax,bathy_map);

% % Ligia min methane composition of 39:50:10.73 percent
% methane:ethane:nitrogen from LeGall 2016
rho_ligiamin = 480.700;
sfct_ligiamin = 2.6E-2;
nu_ligiamin = 6.363E-4;
sigH_ligiamin = makeWaves(windspeeds,rho_ligiamin,sfct_ligiamin,nu_ligiamin,bathy_map);

% sigH_methane = makeWaves(windspeeds,rho_methane,nu_methane,bathy_map);
% sigH_ethane = makeWaves(windspeeds,rho_ethane,nu_ethane,bathy_map);
% sigH_ontario = makeWaves(windspeeds,rho_ontario,nu_ontario,bathy_map);
% sigH_punga = makeWaves(windspeeds,rho_punga,nu_punga,bathy_map);

% figure
% plot(windspeeds,sigH_methane(:,2),'-o')
% hold on
% plot(windspeeds,sigH_ethane(:,2),'-o')
% legend('Methane-N2','Ethane-N2')
% xlabel('wind speed [m/s]')
% ylabel('sig H [cm]')
% hold off

% %% Plotting for varying methane:ethane:nitrogen compositions
% figure
% plot(windspeeds,sigH_0_96_4(:,2),'-o', 'Color', [66/255 0/255 83/255] )
% hold on
% plot(windspeeds,sigH_10_86_4(:,2),'-o', 'Color', [69/255 32/255 118/255])
% plot(windspeeds,sigH_19_76_5(:,2),'-o', 'Color', [62/255 66/255 137/255])
% plot(windspeeds,sigH_28_66_6(:,2),'-o', 'Color', [49/255 96/255 143/255])
% plot(windspeeds,sigH_37_55_8(:,2),'-o', 'Color', [38/255 120/255 144/255])
% plot(windspeeds,sigH_45_45_10(:,2),'-o', 'Color', [26/255 145/255 141/255])
% plot(windspeeds,sigH_53_35_12(:,2),'-o', 'Color', [30/255 169/255 133/255])
% plot(windspeeds,sigH_60_26_14(:,2),'-o', 'Color', [66/255 192/255 114/255])
% plot(windspeeds,sigH_66_16_18(:,2),'-o', 'Color', [123/255 210/255 79/255])
% plot(windspeeds,sigH_71_8_21(:,2),'-o', 'Color', [187/255 224/255 35/255])
% plot(windspeeds,sigH_75_0_25(:,2),'-o', 'Color', [253/255 232/255 33/255])
% legend('0:96:4','10:86:4','19:76:5','28:66:6','37:55:8','45:45:10','53:35:12','60:26:14','66:16:28','71:8:21','75:0:25');
% xlabel('wind speed [m/s]')
% ylabel('sig H [cm]')
% title('windspeed vs sig H at 92K for varying methane:ehtane:nitrogencompositions')
% hold off
% saveas(gcf,'VaryingComps.png')

%% Plotting for varying lake compositions
figure
plot(windspeeds,sigH_ontario(:,2),'-o', 'Color', [66/255 0/255 83/255] )
hold on
plot(windspeeds,sigH_punga(:,2),'-o', 'Color', [38/255 120/255 144/255])
plot(windspeeds,sigH_ligiamax(:,2),'-o', 'Color', [66/255 192/255 114/255])
plot(windspeeds,sigH_ligiamin(:,2),'-o', 'Color', [253/255 232/255 33/255])
legend('Ontario','Punga','Ligia (max)','Ligia (min)');
xlabel('wind speed [m/s]')
ylabel('sig H [cm]')
title('windspeed vs sig H at 92K for varying lake compositions')
hold off
saveas(gcf,'VaryingLakes.png')