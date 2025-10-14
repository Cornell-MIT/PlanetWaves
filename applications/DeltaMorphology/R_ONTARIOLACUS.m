clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Housekeeping

% Ontario Lacus shoreline
addpath(fullfile('..','..','\data\Titan\TitanLakes\shoreline'))
% Juan's winds
addpath(fullfile('..','..','data/Titan/TAMwTopo/')) 

jet_wrap = vertcat(jet,flipud(jet)); % circular colormap
viridis_top = viridis;
viridis_top = viridis_top(1:end/2,:);
plasma_bottom = plasma;
plasma_bottom = plasma_bottom(end/2:end,:);
distinct_cmap = [viridis_top;plasma_bottom];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate max flux from waves

% Load in shoreline
load('ontariolacus_shoreline.mat')
x = lon; y = lati;
x(1) = []; y(1) = [];
[x, y] = make_circle(x,y); % make a simple shape for testing

% Calculate the regional shoreline orientaiton at each point
window_size_ang = 100; % size of window to define the regional shoreline orientation over
shoreline_angle = calc_regional_shoreline_angle(x, y,window_size_ang);


% Specify wave climate
wave_mean = 270; % degrees CCW from E
wave_std = 500;
wave = 0:359;
E_pdf = vonMises(deg2rad(wave_mean),deg2rad(wave_std),deg2rad(wave));

% Calculate wave flux of sediment
Qsmax = wave_flux(E_pdf,rad2deg(shoreline_angle));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTS 
% Plot of shoreline with quiver showing the window size calculating the
% regional angle over
figure('Name','window size for angle estimation');
plot(x,y,'-k','MarkerFaceColor','k')
hold on
quiver(x(1),y(1),x(window_size_ang) - x(1),y(window_size_ang) - y(1),'k')
title('window size')
axis equal padded

% Plot showing the regional shoreline angle at each point
figure('Name','shoreline angle')
scatter(x,y,10,(rad2deg(shoreline_angle)))
colorbar;
title('regional shoreline angle')
axis equal padded

% Plot showing the maximum flux from the waves
figure('Name','Qsmax')
scatter(x,y,10,Qsmax)
colormap(jet)
colorbar
title('Qsmax of waves')
axis equal padded

% 
% 
% Titan_cold.rho_s = 0.95*1000; 
% Titan_cold.rho = 0.67*1000;
% Titan_cold.nu = 0.003/10000;
% Titan_cold.g = 1.35;
% 
% vidflumina.coldbedload_D50 = [3.8 10]/100;
% vidflumina.susload_D50 = 6.35e-5;
% vidflumina.width = [100 175];
% vidflumina.slope = [0.0011 0.0015];
% 
% [river_susload,river_bedload] = riverine_flux(Titan_cold.rho_s,Titan_cold.rho,Titan_cold.nu,Titan_cold.g,vidflumina.coldbedload_D50(1),vidflumina.susload_D50,vidflumina.width(1),NaN);





ATM_2_PASCAL = 101325;
% RRR = 8.314; 
% 
% planet.rho_liquid = 588.15;                                            % TITANPOOL   
% planet.nu_liquid = 8.084e-7;                                           % TITANPOOL   
% planet.nua = 6.4e-9;                                                   % vapor N2 viscocity at 92 K (DIPPR)
% planet.gravity = 1.352;                                                % Titan gravity
% planet.surface_temp = 92;                                              % ~avg Titan surface temp
% planet.surface_press = 1.5*ATM_2_PASCAL;                               % ~avg Titan surface pressure
% planet.surface_tension = 0.032766;                                     % TITANPOOL
% planet.kgmolwt = 0.028;                                                % (N2)
% 
% rho_a = planet.surface_press*planet.kgmolwt/(RRR*planet.surface_temp);
% U10 = 0.4:0.1:4;
% X = [0:10:180].*1000;
% for i = 1:numel(U10)
%     for j = 1:numel(X)
%         H_sig(i,j) = (planet.nu_liquid^(0.002)/(planet.nua^(0.003)))*((rho_a/planet.rho_liquid)^0.36)*(((X(j)^0.07)*(U10(i)^1.6))/(planet.gravity^0.7));
%     end
% end
% 
% figure('Name','sig wave height from eqn')
% contourf(X./1000,U10,H_sig);
% colorbar
% ylabel('wind speed m/s')
% xlabel('Fetch km')



