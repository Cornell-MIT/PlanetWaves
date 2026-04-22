clc
clear
close all

%
addpath(genpath(fullfile('subroutines')));
run_waves = 0; % if true, will run PlanetWaves (Titan_DeepwaterWaves.m), otherwise will load wave results from past run

% wave model
addpath(genpath(fullfile('..','..','planetwaves')))
% GCM surface winds
addpath(fullfile('..','..','data/Titan/TAMwTopo/'))
% simple bathymetry for Ontario Lacus using constant slope
addpath(fullfile('..','..','data/Titan/TitanLakes/Bathymetries/bathtub_bathy'))


if run_waves
    Titan_DeepwaterWaves;
else
    pathname = [pwd,'\past_runs\Titan_DeepWaterWaves.mat'];
    last_created = dir(pathname).date;
    fprintf('Loading data from file run: %s\n',last_created);
    % Results from last run of Titan_DeepwaterWaves
    load('./past_runs/Titan_DeepwaterWaves.mat','test_speeds','time_to_run','wind_direction','zDep','buoy_loc','lakes','Planet','Model','H0','L0','T','C0','Cg0')
end

lakecolors = {'#F2D1C9','#E086D3','#8332AC','#462749'};

figure('Name','u v H inputs');
for i = 1:numel(lakecolors)
    plot(test_speeds,H0(i,:),'-','LineWidth',5,'Color',lakecolors{i},'DisplayName',lakes{i})
    hold on;
end
xlabel('wind speed [m/s]')
ylabel('wave height [m]')
grid on
legend('show')


% TAM
load('OL_winds.mat','mag_wind','angle_wind')
wind_mag = mag_wind;
wind_angle_deg = wrapTo360(round((angle_wind))); % wind going from west to east = 0, positive CCW

% figure;
% wind_rose(wind_angle_deg,wind_mag)
% title('TAM, Ontario lacus')

H_wave = H0(2,:);
L_wave = L0(2,:);
T_wave = T(2,:);

g = 1.352;
rho_s = 940;
rho = 588.15; 
R = rho_s/rho - 1;
nu = 8.084e-7;
beta = 1e-3;
alpha0 = wrapToPi(deg2rad(0:359));
D50 = 6.35e-5; % d50 = 6.35e-5:1e-5:0.1;  % [fines gravel]

linecolor = flip(winter(numel(test_speeds)));
figure;
hold on;
for i = 1:numel(test_speeds)
    Qs(i,:) = calc_Diegaard_wave_flux(H_wave(i), L_wave(i), rho_s,rho, g, D50, nu, beta, alpha0);
    Qs_CERC(i,:) = CERC(H_wave(i),T_wave(i),rho,rho_s,g,rad2deg(alpha0));
    plot(rad2deg(alpha0),Qs(i,:),'ob','LineWidth',3,'color',linecolor(i,:))
    %plot(rad2deg(alpha0),Qs_CERC(i,:),'.r')
    [Qs_max(i),i_angle] = max(Qs(i,:));
    %plot(rad2deg(alpha0(i_angle)),Qs_max(i),'or','MarkerFaceColor','r')
end
colormap(flip(winter(i)));
cb = colorbar;
caxis([0 max(test_speeds)])
cb.Label.String = 'u m/s';
grid on;
title('Diegaard')

xlabel('wave approach angle $\alpha_0$ [deg]','Interpreter','latex')
ylabel('sediment flux Q$_s$ [m$^3$/s]','Interpreter','latex')


bin_size = 0.1;
figHandle = figure;
yyaxis left
h = histogram(wind_mag,'BinEdges',0:bin_size:max(wind_mag)+bin_size,'Normalization','probability','DisplayStyle','stairs','EdgeColor','b','LineWidth',2);
axisHandle = figHandle.Children;
histHandle = axisHandle.Children;
histHandle.BinEdges = histHandle.BinEdges + histHandle.BinWidth/2;
aa = histHandle.BinEdges;
percent_time = h.Values;
percent_time = [0 percent_time];
u_speed = aa - bin_size/2;
grid on
xlabel('wind speed [m/s]')
ylabel('PDF')
yyaxis right
Qs_normalized = Qs_max./max(Qs_max);
Qs_normalized(isnan(Qs_normalized)) = 0;
Qs_normalized = [0 0 0 Qs_normalized]; % no Qs for 0,0.1,0.2 m/s 
plot(u_speed,Qs_normalized,'--','Color','r','LineWidth',2)
hold on;
Wolman_Miller = Qs_normalized.*percent_time;
%Wolman_Miller = ones(size(percent_time)).*percent_time;
Wolman_Miller = Wolman_Miller./max(Wolman_Miller);
plot(u_speed,Wolman_Miller,'-','LineWidth',5,'Color','m')
%plot(test_speeds,H0(i,:),'-','LineWidth',5,'Color',lakecolors{i},'DisplayName',lakes{i})
ylabel('Q_s/Q_{s,max}')
ax = gca;
xlim([0 5])
ax.FontSize = 20; 

[~,I] = max(Wolman_Miller);
most_effect_wind = u_speed(I);

fprintf('The most geomorpholigically effective wind speed is %0.2f\n',most_effect_wind)
% Calculate reoccurence interval
[r,x_sorted] = calc_recurrence_interval(mag_wind,most_effect_wind);