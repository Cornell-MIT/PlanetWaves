clc
clear
close all

% make plots of signifigant wave height for different compositions

m = 31;                                                                    % number of gridpoints in x-direction
n = 15;                                                                    % number of gridpoints in y-direction

windspeeds = [0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5]; %for time_step_size = 100 and num_time_steps = 200, time per run = 50 min, total time = 14 hr 30 min
%windspeeds = [1 2 3]; %for time_step_size = 100 and num_time_steps = 200, time per run = 5 min, total time = 1 hr 20 min
wind_dir = 0;
planet_gravity = 1.35;
planet_temp = 92;
planet_press = 1.5*101300;
bathy_map = 100.*ones(m,n);
gridX = m;
gridY = n;
time_step_size = 100; %usually is 100
num_time_steps = 200; %usually is 200


%% methane:ethane:nitrogen from titanpool at 92K (steckloff et al., 2020) and NIST
% 0:96:4.27
rho_0_96_4 = 652.735;
sfct_0_96_4 = 1.209E-3;
nu_0_96_4 = 1.209E-3/rho_0_96_4;
sigH_0_96_4 = makeWaves(windspeeds,wind_dir,rho_0_96_4,nu_0_96_4, ...
    planet_gravity,planet_temp,planet_press,sfct_0_96_4,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_0_96_4(size(sigH_0_96_4,1),:), 'LineWidth',2 ) ;
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for 0:96:4.27 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'0_96_4_timevSigH.png')
figure
sigH_0_96_4_wind = sigH_0_96_4(:,~any(isnan(sigH_0_96_4)));
plot(windspeeds,sigH_0_96_4_wind(:,size(sigH_0_96_4_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for 0:96:4.27 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'0_96_4_windvSigH.png')

% 10:86:4.96
rho_10_86_4 = 641.081;
sfct_10_86_4 = 3.562E-2;
nu_10_86_4 = 1.106E-3/rho_10_86_4;
sigH_10_86_4 = makeWaves(windspeeds,wind_dir,rho_10_86_4,nu_10_86_4, ...
    planet_gravity,planet_temp,planet_press,sfct_10_86_4,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_10_86_4(size(sigH_10_86_4,1),:), 'LineWidth',2  );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for 10:86:4.96 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'0_96_4_timevSigH.png')
figure
sigH_10_86_4_wind = sigH_10_86_4(:,~any(isnan(sigH_10_86_4)));
plot(windspeeds,sigH_10_86_4_wind(:,size(sigH_10_86_4_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for 10:86:4.96 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'10_86_4_windvSigH.png')

% 19:76:5.53
rho_19_76_5 = 626.643;
sfct_19_76_5 = 3.353E-2;
nu_19_76_5 = 1.004E-3/rho_19_76_5;
sigH_19_76_5 = makeWaves(windspeeds,wind_dir,rho_19_76_5,nu_19_76_5, ...
    planet_gravity,planet_temp,planet_press,sfct_19_76_5,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_19_76_5(size(sigH_19_76_5,1),:), 'LineWidth',2  );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for 19:76:5.53 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'19_76_5_timevSigH.png')
figure
sigH_19_76_5_wind = sigH_19_76_5(:,~any(isnan(sigH_19_76_5)));
plot(windspeeds,sigH_19_76_5_wind(:,size(sigH_19_76_5_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for 19:76:5.53 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'19_76_5_windvSigH.png')


% 28:66:6.35
rho_28_66_6 = 612.61;
sfct_28_66_6 = 3.144E-2;
nu_28_66_6 = 9.02E-4/rho_28_66_6;
sigH_28_66_6 = makeWaves(windspeeds,wind_dir,rho_28_66_6,nu_28_66_6, ...
    planet_gravity,planet_temp,planet_press,sfct_28_66_6,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_28_66_6(size(sigH_28_66_6,1),:), 'LineWidth',2  );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for 28:66:6.35 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'28_66_6_timevSigH.png')
figure
sigH_28_66_6_wind = sigH_28_66_6(:,~any(isnan(sigH_28_66_6)));
plot(windspeeds,sigH_28_66_6_wind(:,size(sigH_28_66_6_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for 28:66:6.35 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'28_66_6_windvSigH.png')

% 37:55:7.61
rho_37_55_8 = 597.556;
sfct_37_55_8 = 2.935E-2;
nu_37_55_8 = 7.998E-4/rho_37_55_8;
sigH_37_55_8 = makeWaves(windspeeds,wind_dir,rho_37_55_8,nu_37_55_8, ...
    planet_gravity,planet_temp,planet_press,sfct_37_55_8,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_37_55_8(size(sigH_37_55_8,1),:), 'LineWidth',2  );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for 37:55:7.61 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'37_55_8_timevSigH.png')
figure
sigH_37_55_8_wind = sigH_37_55_8(:,~any(isnan(sigH_37_55_8)));
plot(windspeeds,sigH_37_55_8_wind(:,size(sigH_37_55_8_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for 37:55:7.61 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'37_55_8_windvSigH.png')

% 45:45:9.42
rho_45_45_10 = 582.749;
sfct_45_45_10 = 2.726E-2;
nu_45_45_10 = 6.976E-4/rho_45_45_10;
sigH_45_45_10 = makeWaves(windspeeds,wind_dir,rho_45_45_10,nu_45_45_10, ...
    planet_gravity,planet_temp,planet_press,sfct_45_45_10,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_45_45_10(size(sigH_45_45_10,1),:), 'LineWidth',2  );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for 45:45:9.42 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'45_45_10_timevSigH.png')
figure
sigH_45_45_10_wind = sigH_45_45_10(:,~any(isnan(sigH_45_45_10)));
plot(windspeeds,sigH_45_45_10_wind(:,size(sigH_45_45_10_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for 45:45:9.42 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'45_45_10_windvSigH.png')

% 53:35:11.52
rho_53_35_12 = 567.498;
sfct_53_35_12 = 2.517E-2;
nu_53_35_12 = 5.954E-4/rho_53_35_12;
sigH_53_35_12 = makeWaves(windspeeds,wind_dir,rho_53_35_12,nu_53_35_12, ...
    planet_gravity,planet_temp,planet_press,sfct_53_35_12,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_53_35_12(size(sigH_53_35_12,1),:), 'LineWidth',2  );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for 53:35:11.52 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'53_35_12_timevSigH.png')
figure
sigH_53_35_12_wind = sigH_53_35_12(:,~any(isnan(sigH_53_35_12)));
plot(windspeeds,sigH_53_35_12_wind(:,size(sigH_53_35_12_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for 53:35:11.52 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'53_35_12_windvSigH.png')

% 60:26:14.27
rho_60_26_14 = 552.24;
sfct_60_26_14 = 2.308E-2;
nu_60_26_14 = 4.932E-4/rho_60_26_14;
sigH_60_26_14 = makeWaves(windspeeds,wind_dir,rho_60_26_14,nu_60_26_14, ...
    planet_gravity,planet_temp,planet_press,sfct_60_26_14,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_60_26_14(size(sigH_60_26_14,1),:), 'LineWidth',2  );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for 60:26:14.27 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'60_26_14_timevSigH.png')
figure
sigH_60_26_14_wind = sigH_60_26_14(:,~any(isnan(sigH_60_26_14)));
plot(windspeeds,sigH_60_26_14_wind(:,size(sigH_60_26_14_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for 60:26:14.27 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'60_26_14_windvSigH.png')

% 66:16:17.70 
rho_66_16_18 = 538.596;
sfct_66_16_18 = 2.098E-2;
nu_66_16_18 = 3.910E-4/rho_66_16_18;
sigH_66_16_18 = makeWaves(windspeeds,wind_dir,rho_66_16_18,nu_66_16_18, ...
    planet_gravity,planet_temp,planet_press,sfct_66_16_18,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_66_16_18(size(sigH_66_16_18,1),:), 'LineWidth',2  );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for 66:16:17.70 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'66_16_18_timevSigH.png')
figure
sigH_66_16_18_wind = sigH_66_16_18(:,~any(isnan(sigH_66_16_18)));
plot(windspeeds,sigH_66_16_18_wind(:,size(sigH_66_16_18_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for 66:16:17.70 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'66_16_18_windvSigH.png')

% 71:8:21.59
rho_71_8_21 = 524.596;
sfct_71_8_21 = 1.889E-2;
nu_71_8_21 = 1.106E-3/rho_71_8_21;
sigH_71_8_21 = makeWaves(windspeeds,wind_dir,rho_71_8_21,nu_71_8_21, ...
    planet_gravity,planet_temp,planet_press,sfct_71_8_21,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_71_8_21(size(sigH_71_8_21,1),:), 'LineWidth',2 );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for 71_8_21 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'71_8_21_timevSigH.png')
figure
sigH_71_8_21_wind = sigH_71_8_21(:,~any(isnan(sigH_71_8_21)));
plot(windspeeds,sigH_71_8_21_wind(:,size(sigH_71_8_21_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for 71:8:21 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'71_8_21_windvSigH.png')

% 75:0:25.24
rho_75_0_25 = 512.577;
sfct_75_0_25 = 1.68E-2;
nu_75_0_25 = 1.867E-4/rho_75_0_25;
sigH_75_0_25 = makeWaves(windspeeds,wind_dir,rho_75_0_25,nu_75_0_25, ...
    planet_gravity,planet_temp,planet_press,sfct_75_0_25,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_75_0_25(size(sigH_75_0_25,1),:), 'LineWidth',2 );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for 75:0:25.24 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'75_0_25_timevSigH.png')
figure
sigH_75_0_25_wind = sigH_75_0_25(:,~any(isnan(sigH_75_0_25)));
plot(windspeeds,sigH_75_0_25_wind(:,size(sigH_75_0_25_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for 75:0:25.24 (methane:ethane:nitrogen)','FontSize',25,'FontWeight','bold')
saveas(gcf,'75_0_25_windvSigH.png')

%% lake compositions from NIST and sources (80 min for 2 windspeeds; 160 min for 4 windspeeds)
% % Ontario composition of 51:38:11 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
rho_ontario = 476.962;
sfct_ontario = 2.579E-2;
nu_ontario = 6.261E-4/rho_ontario;
sigH_ontario = makeWaves(windspeeds,wind_dir,rho_ontario,nu_ontario, ...
    planet_gravity,planet_temp,planet_press,sfct_ontario,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps); 
figure
plot(1:num_time_steps,sigH_ontario(size(sigH_ontario,1),:), 'LineWidth',2  );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for Ontario','FontSize',25,'FontWeight','bold')
saveas(gcf,'Ontari_timevSigH.png')
figure
sigH_ontario_wind = sigH_ontario(:,~any(isnan(sigH_ontario)));
plot(windspeeds,sigH_ontario_wind(:,size(sigH_ontario_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for Ontario','FontSize',25,'FontWeight','bold')
saveas(gcf,'Ontario_windvSigH.png')

% % Punga composition of 80:0:20 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
rho_punga = 361.006;
sfct_punga = 1.680E-2;
nu_punga = 1.866E-4/rho_punga;
sigH_punga = makeWaves(windspeeds,wind_dir,rho_punga,nu_punga, ...
    planet_gravity,planet_temp,planet_press,sfct_punga,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_punga(size(sigH_punga,1),:), 'LineWidth',2 );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for Punga','FontSize',25,'FontWeight','bold')
saveas(gcf,'Punga_timevSigH.png')
figure
sigH_punga_wind = sigH_punga(:,~any(isnan(sigH_punga)));
plot(windspeeds,sigH_punga_wind(:,size(sigH_punga_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for Punga','FontSize',25,'FontWeight','bold')
saveas(gcf,'Punga_windvSigH.png')

% % Ligia max methane composition of 100:0:0 percent
% methane:ethane:nitrogen from LeGall 2016
rho_ligiamax = 449.83;
sfct_ligiamax = 1.68E-2;
nu_ligiamax = 1.866E-4/rho_ligiamax;
sigH_ligiamax = makeWaves(windspeeds,wind_dir,rho_ligiamax,nu_ligiamax, ...
    planet_gravity,planet_temp,planet_press,sfct_ligiamax,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_ligiamax(size(sigH_ligiamax,1),:), 'LineWidth',2 );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for Ligia (max methane)','FontSize',25,'FontWeight','bold')
saveas(gcf,'LigiaMax_timevSigH.png')
figure
sigH_ligiamax_wind = sigH_ligiamax(:,~any(isnan(sigH_ligiamax)));
plot(windspeeds,sigH_ligiamax_wind(:,size(sigH_ligiamax_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for Ligia (max methane)','FontSize',25,'FontWeight','bold')
saveas(gcf,'LigiaMax_windvSigH.png')

% % Ligia min methane composition of 39:50:10.73 percent
% methane:ethane:nitrogen from LeGall 2016
rho_ligiamin = 480.700;
sfct_ligiamin = 2.6E-2;
nu_ligiamin = 6.363E-4/rho_ligiamin;
sigH_ligiamin = makeWaves(windspeeds,wind_dir,rho_ligiamin,nu_ligiamin, ...
    planet_gravity,planet_temp,planet_press,sfct_ligiamin,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps); 
figure
plot(1:num_time_steps,sigH_ligiamin(size(sigH_ligiamin,1),:), 'LineWidth',2 );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for Ligia (min methane)','FontSize',25,'FontWeight','bold')
saveas(gcf,'LigiaMin_timevSigH.png')
figure
sigH_ligiamin_wind = sigH_ligiamin(:,~any(isnan(sigH_ligiamin)));
plot(windspeeds,sigH_ligiamin_wind(:,size(sigH_ligiamin_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for Ligia (min methane)','FontSize',25,'FontWeight','bold')
saveas(gcf,'LigiaMin_windvSigH.png')

% % Kraken methane composition of 70:16:14 percent
% methane:ethane:nitrogen from Poggiali et al., 2020
rho_kraken = 527.625;
sfct_kraken = 2.07E-2;
nu_kraken = 3.81E-4/rho_kraken;
sigH_kraken = makeWaves(windspeeds,wind_dir,rho_kraken,nu_kraken, ...
    planet_gravity,planet_temp,planet_press,sfct_kraken,bathy_map, ...
    gridX,gridY,time_step_size,num_time_steps);
figure
plot(1:num_time_steps,sigH_kraken(size(sigH_kraken,1),:), 'LineWidth',2 );
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for Kraken','FontSize',25,'FontWeight','bold')
saveas(gcf,'Kraken_timevSigH.png') 
figure
sigH_kraken_wind = sigH_kraken(:,~any(isnan(sigH_kraken)));
plot(windspeeds,sigH_kraken_wind(:,size(sigH_kraken_wind,2)), 'LineWidth',2 );
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for Kraken','FontSize',25,'FontWeight','bold')
saveas(gcf,'Kraken_windvSigH.png')

%% Plotting for varying lake compositions
figure
plot(1:num_time_steps,sigH_0_96_4(size(sigH_0_96_4,1),:),'-o', 'Color', [66/255 0/255 83/255],'LineWidth',2  )
hold on
plot(1:num_time_steps,sigH_10_86_4(size(sigH_10_86_4,1),:),'-o', 'Color', [69/255 32/255 118/255],'LineWidth',2 )
plot(1:num_time_steps,sigH_19_76_5(size(sigH_19_76_5,1),:),'-o', 'Color', [62/255 66/255 137/255],'LineWidth',2 )
plot(1:num_time_steps,sigH_28_66_6(size(sigH_28_66_6,1),:),'-o', 'Color', [49/255 96/255 143/255],'LineWidth',2 )
plot(1:num_time_steps,sigH_37_55_8(size(sigH_37_55_8,1),:),'-o', 'Color', [38/255 120/255 144/255],'LineWidth',2 )
plot(1:num_time_steps,sigH_45_45_10(size(sigH_45_45_10,1),:),'-o', 'Color', [26/255 145/255 141/255],'LineWidth',2 )
plot(1:num_time_steps,sigH_53_35_12(size(sigH_53_35_12,1),:),'-o', 'Color', [30/255 169/255 133/255],'LineWidth',2 )
plot(1:num_time_steps,sigH_60_26_14(size(sigH_60_26_14,1),:),'-o', 'Color', [66/255 192/255 114/255],'LineWidth',2 )
plot(1:num_time_steps,sigH_66_16_18(size(sigH_66_16_18,1),:),'-o', 'Color', [123/255 210/255 79/255],'LineWidth',2 )
plot(1:num_time_steps,sigH_71_8_21(size(sigH_71_8_21,1),:),'-o', 'Color', [187/255 224/255 35/255],'LineWidth',2 )
plot(1:num_time_steps,sigH_75_0_25(size(sigH_75_0_25,1),:),'-o', 'Color', [253/255 232/255 33/255],'LineWidth',2 )
legend('0:96:4','10:86:4','19:76:5','28:66:6','37:55:8','45:45:10','53:35:12','60:26:14','66:16:28','71:8:21','75:0:25','Location','northwest','FontSize',12);
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for varying methane:ehtane:nitrogen compositions','FontSize',25,'FontWeight','bold')
hold off
saveas(gcf,'VaryingComps_timevSigH.png')

figure
plot(windspeeds,sigH_0_96_4_wind(:,size(sigH_0_96_4_wind,2)),'-o', 'Color', [66/255 0/255 83/255] )
hold on
plot(windspeeds,sigH_10_86_4_wind(:,size(sigH_10_86_4_wind,2)),'-o', 'Color', [69/255 32/255 118/255],'LineWidth',2 )
plot(windspeeds,sigH_19_76_5_wind(:,size(sigH_19_76_5_wind,2)),'-o', 'Color', [62/255 66/255 137/255],'LineWidth',2 )
plot(windspeeds,sigH_28_66_6_wind(:,size(sigH_28_66_6_wind,2)),'-o', 'Color', [49/255 96/255 143/255],'LineWidth',2 )
plot(windspeeds,sigH_37_55_8_wind(:,size(sigH_37_55_8_wind,2)),'-o', 'Color', [38/255 120/255 144/255],'LineWidth',2 )
plot(windspeeds,sigH_45_45_10_wind(:,size(sigH_45_45_10_wind,2)),'-o', 'Color', [26/255 145/255 141/255],'LineWidth',2 )
plot(windspeeds,sigH_53_35_12_wind(:,size(sigH_53_35_12_wind,2)),'-o', 'Color', [30/255 169/255 133/255],'LineWidth',2 )
plot(windspeeds,sigH_60_26_14_wind(:,size(sigH_60_26_14_wind,2)),'-o', 'Color', [66/255 192/255 114/255],'LineWidth',2 )
plot(windspeeds,sigH_66_16_18_wind(:,size(sigH_66_16_18_wind,2)),'-o', 'Color', [123/255 210/255 79/255],'LineWidth',2 )
plot(windspeeds,sigH_71_8_21_wind(:,size(sigH_71_8_21_wind,2)),'-o', 'Color', [187/255 224/255 35/255],'LineWidth',2 )
plot(windspeeds,sigH_75_0_25_wind(:,size(sigH_75_0_25_wind,2)),'-o', 'Color', [253/255 232/255 33/255],'LineWidth',2 )
legend('0:96:4','10:86:4','19:76:5','28:66:6','37:55:8','45:45:10','53:35:12','60:26:14','66:16:28','71:8:21','75:0:25','Location','northwest','FontSize',12);
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for methane:ehtane:nitrogen compositions','FontSize',25,'FontWeight','bold')
hold off
saveas(gcf,'VaryingComps_windvSigH.png')

%% Plotting for different lakes
figure
plot(1:num_time_steps,sigH_ontario(size(sigH_ontario,1),:),'-o', 'Color', [66/255 0/255 83/255],'LineWidth',2 )
hold on
plot(1:num_time_steps,sigH_punga(size(sigH_punga,1),:),'-o', 'Color', [62/255 66/255 137/255],'LineWidth',2 );
plot(1:num_time_steps,sigH_ligiamax(size(sigH_ligiamax,1),:),'-o', 'Color', [38/255 120/255 144/255],'LineWidth',2 );
plot(1:num_time_steps,sigH_ligiamin(size(sigH_ligiamin,1),:),'-o', 'Color', [66/255 192/255 114/255],'LineWidth',2 );
plot(1:num_time_steps,sigH_kraken(size(sigH_kraken,1),:),'-o', 'Color', [253/255 232/255 33/255],'LineWidth',2 );
legend('Ontario','Punga','Ligia (max)','Ligia (min)','Kraken','Location','northwest','FontSize',12);
xlabel('# time steps (delt = .1 s)','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('timestep vs sig H at 92K for varying lake compositions','FontSize',25,'FontWeight','northwest')
hold off
saveas(gcf,'TitanLakes_timevSigH.png')

figure
plot(windspeeds,sigH_ontario_wind(:,size(sigH_ontario_wind,2)),'-o', 'Color', [66/255 0/255 83/255],'LineWidth',2 )
hold on
plot(windspeeds,sigH_punga_wind(:,size(sigH_punga_wind,2)),'-o', 'Color', [62/255 66/255 137/255],'LineWidth',2 );
plot(windspeeds,sigH_ligiamax_wind(:,size(sigH_ligiamax_wind,2)),'-o', 'Color', [38/255 120/255 144/255],'LineWidth',2 );
plot(windspeeds,sigH_ligiamin_wind(:,size(sigH_ligiamin_wind,2)),'-o', 'Color', [66/255 192/255 114/255],'LineWidth',2 );
plot(windspeeds,sigH_kraken_wind(:,size(sigH_kraken_wind,2)),'-o', 'Color', [253/255 232/255 33/255],'LineWidth',2 );
legend('Ontario','Punga','Ligia (max)','Ligia (min)','Kraken','Location','northwest','FontSize',12);
xlabel('windspeed [m/s]','FontSize',16,'FontWeight','bold')
ylabel('sig H [m]','FontSize',16,'FontWeight','bold')
xlim([1 size(sigH_0_96_4,1)]);
ylim([0 1]);
title('windspeed vs sig H at 92K for varying lakes','FontSize',25,'FontWeight','bold')
hold off
saveas(gcf,'TitanLakes_windvSigH.png')