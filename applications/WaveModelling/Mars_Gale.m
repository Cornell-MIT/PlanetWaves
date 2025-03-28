clc
clear
close all

addpath(fullfile('..','..','data','Mars','StevensRubin2022'))
addpath(fullfile('..','..','planetwaves'))
addpath(fullfile('..','..','planetwaves','pre_analysis/'))

max_depth = 10;

%downwind
x = linspace(0, 20000, 50); 
% crosswind
y = 1:5; 

depth = zeros(length(y), length(x));

% slope is discontinous
fetch = 20*1000;
slope_break = 18800;
flat_region = x <= slope_break; 
slope_region = x > slope_break;

% depth only varying in downwind
depth(:, flat_region) = max_depth; 

% Linear slope in the last 1.2 km 
depth(:, slope_region) = repmat(max_depth * (1 - (x(slope_region) - slope_break) / (fetch-slope_break)), length(y), 1);

% Create meshgrid for plotting
[X, Y] = meshgrid(x / 1000, y); % Convert to km

% Plot surface
figure;
h = pcolor(X, Y, depth);
set(h,'EdgeColor','none');
xlabel('downwind (km)');
ylabel('crosswind (km)');
zlabel('Depth (m)');
title('Gale (Rubin+2022)');
colorbar;
hold on;
xline(0:1.2:18)


planet_to_run = 'Mars-high';
time_to_run = 2*60;
test_speeds = [1.8 2.15 4.3 8.6 17.2];
wind_direction = 0;
buoy_loc = [48 2];

[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,depth,buoy_loc);
Model.gridX = 400;                                              
Model.gridY = 10*1000;  
Model.Fdim = 24;
Model.min_freq = 0.05;                                                     % minimum frequency to model
Model.max_freq = 1;  
Planet.surface_press = 60000; % 600 mbar
Planet.surface_temp = 273;
make_input_map(Planet,Model,Wind)

for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);
    Model = calc_cutoff_freq(Planet,Model,Wind);

    [myHsig{i}, htgrid{i}, wn_e_spectrum, ~ , ~ , ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
    if ~isempty(wn_e_spectrum{end})
        energy{i} = squeeze(sum(wn_e_spectrum{end}.E(Model.long,Model.lat,:,:),4));
        wn{i} = squeeze(sum(wn_e_spectrum{end}.k(Model.long,Model.lat,:,:),4));
        cg{i} = squeeze(sum(wn_e_spectrum{end}.cg(Model.long,Model.lat,:,:),4));
        if numel(htgrid{i}{end}) > 0
            h_v_fetch(i,:) = htgrid{i}{end}(:,buoy_loc(2));
        else
            h_v_fetch(i,:) = NaN(50,1);
        end
    end
end
%make_plots(Planet,Model,Wind,test_speeds,myHsig, htgrid,energy,wn)


load('stevensrubins2022.mat')
x_SWAN = (marswaveparams.distance_km);
H_SWAN = (marswaveparams.hsig_m);

% p6_u7p5 = marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==7.5;
% x_SWAN_p6_u7p5 = x_SWAN(p6_u7p5);
% H_SWAN_p6_u7p5 = H_SWAN(p6_u7p5);
% 
% p6_u15 = marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==15;
% x_SWAN_p6_u15 = x_SWAN(p6_u15);
% H_SWAN_p6_u15 = H_SWAN(p6_u15);
% 
% p6_u30 = marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==30;
% x_SWAN_p6_u30 = x_SWAN(p6_u30);
% H_SWAN_p6_u30 = H_SWAN(p6_u30);
% 
% p6_u60 = marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==60;
% x_SWAN_p6_u60 = x_SWAN(p6_u60);
% H_SWAN_p6_u60 = H_SWAN(p6_u60);
% 
% p6_u120 = marswaveparams.pressure_mbar==6 & marswaveparams.wind_speed_ms==120;
% x_SWAN_p6_u120 = x_SWAN(p6_u120);
% H_SWAN_p6_u120 = H_SWAN(p6_u120);
% 
% figure;
% plot(x_SWAN_p6_u7p5,H_SWAN_p6_u7p5);
% hold on
% plot(x_SWAN_p6_u15,H_SWAN_p6_u15);
% plot(x_SWAN_p6_u30,H_SWAN_p6_u30);
% plot(x_SWAN_p6_u60,H_SWAN_p6_u60);
% plot(x_SWAN_p6_u120,x_SWAN_p6_u120);
% 
% 
% p60_u3p67 = marswaveparams.pressure_mbar==60 & marswaveparams.wind_speed_ms==3.67;
% x_SWAN_p60_u3p67 = x_SWAN(p60_u3p67);
% H_SWAN_p60_u3p67 = H_SWAN(p60_u3p67);
% 
% p60_u7p35 = marswaveparams.pressure_mbar==60 & marswaveparams.wind_speed_ms==7.35;
% x_SWAN_p60_u7p35 = x_SWAN(p60_u7p35);
% H_SWAN_p60_u7p35 = H_SWAN(p60_u7p35);
% 
% p60_u14p7 = marswaveparams.pressure_mbar==60 & marswaveparams.wind_speed_ms==14.7;
% x_SWAN_p60_u14p7 = x_SWAN(p60_u14p7);
% H_SWAN_p60_u14p7 = H_SWAN(p60_u14p7);
% 
% p60_u29p4 = marswaveparams.pressure_mbar==60 & marswaveparams.wind_speed_ms==29.4;
% x_SWAN_p60_u29p4 = x_SWAN(p60_u29p4);
% H_SWAN_p60_u29p4 = H_SWAN(p60_u29p4);
% 
% p60_u58p8 = marswaveparams.pressure_mbar==60 & marswaveparams.wind_speed_ms==58.8;
% x_SWAN_p60_u58p8 = x_SWAN(p60_u58p8);
% H_SWAN_p60_u58p8 = H_SWAN(p60_u58p8);
% 
% 
% plot(x_SWAN_p60_u3p67,H_SWAN_p60_u3p67,'--');
% plot(x_SWAN_p60_u7p35,H_SWAN_p60_u7p35,'--');
% plot(x_SWAN_p60_u14p7,H_SWAN_p60_u14p7,'--');
% plot(x_SWAN_p60_u29p4,H_SWAN_p60_u29p4,'--');
% plot(x_SWAN_p60_u58p8,H_SWAN_p60_u58p8,'--');

clr = ['r','g','b','c','m'];
figure;
hold on;
for i = 1:numel(test_speeds)
    plot(x./1000,h_v_fetch(i,:),'-','Color',clr(i),'LineWidth',3)
 
end

p600_u1p08 = marswaveparams.pressure_mbar==600 & marswaveparams.wind_speed_ms==1.08;
x_SWAN_p600_u1p08 = x_SWAN(p600_u1p08);
H_SWAN_p600_u1p08 = H_SWAN(p600_u1p08);

p600_u2p15 = marswaveparams.pressure_mbar==600 & marswaveparams.wind_speed_ms==2.15;
x_SWAN_p600_u2p15 = x_SWAN(p600_u2p15);
H_SWAN_p600_u2p15 = H_SWAN(p600_u2p15);

p600_u4p3 = marswaveparams.pressure_mbar==600 & marswaveparams.wind_speed_ms==4.3;
x_SWAN_p600_u4p3 = x_SWAN(p600_u4p3);
H_SWAN_p600_u4p3 = H_SWAN(p600_u4p3);

p600_u8p6 = marswaveparams.pressure_mbar==600 & marswaveparams.wind_speed_ms==8.6;
x_SWAN_p600_u8p6 = x_SWAN(p600_u8p6);
H_SWAN_p600_u8p6 = H_SWAN(p600_u8p6);

p600_u17p2 = marswaveparams.pressure_mbar==600 & marswaveparams.wind_speed_ms==17.2;
x_SWAN_p600_u17p2 = x_SWAN(p600_u17p2);
H_SWAN_p600_u17p2 = H_SWAN(p600_u17p2);

plot(x_SWAN_p600_u1p08,H_SWAN_p600_u1p08,':','Color',clr(1),'LineWidth',3);
plot(x_SWAN_p600_u2p15,H_SWAN_p600_u2p15,':','Color',clr(2),'LineWidth',3);
plot(x_SWAN_p600_u4p3,H_SWAN_p600_u4p3,':','Color',clr(3),'LineWidth',3);
plot(x_SWAN_p600_u8p6,H_SWAN_p600_u8p6,':','Color',clr(4),'LineWidth',3);
plot(x_SWAN_p600_u17p2,H_SWAN_p600_u17p2,':','Color',clr(5),'LineWidth',3);
title('P = 600 mbar')
xlabel('fetch [km]')
ylabel('wave height [m]')



% legend('u = 7.5','u = 15','u = 30','u = 60','u = 120','u = 3.67','u = 7.35','u = 14.7','u = 29.4','u = 58.8','u = 1.08','u = 2.15','u = 4.3','u = 8.6','u = 17.2')



