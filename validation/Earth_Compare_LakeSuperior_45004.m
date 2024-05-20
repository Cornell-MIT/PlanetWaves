clc
clear
close all

% MODELED VS OBSERVED AT LAKE SUPERIOR DEEPWATER BUOY

addpath(fullfile('..','planetwaves'))  
addpath(fullfile('..','planetwaves','pre_analysis'))
addpath(fullfile('..','data','Earth'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT OBSERVATION 
loc = fullfile('..','data','Earth','WindFetchLS_45004.csv');
A = readtable(loc,'VariableNamingRule', 'preserve');
A.Properties.VariableNames = {'Dir_deg','Fetch_m'};
dir_to_fetch = containers.Map(A.Dir_deg,A.Fetch_m);


[u,h,angle] = table_quiet_times();

WIND_SPEED = [];
WIND_FETCH = [];
WIND_HEIGHT = [];

for i = 1:numel(u)
    
    if numel(u{i}) > 1
        WIND_SPEED = [WIND_SPEED; u{i}];
        WIND_HEIGHT = [WIND_HEIGHT; h{i}];
        for j = 1:numel(u{i})
            WIND_FETCH = [WIND_FETCH; dir_to_fetch(angle{i}(j))];
        end
    end

end


Events = table(WIND_SPEED,WIND_FETCH,WIND_HEIGHT);
Events.Properties.VariableNames = {'U','F','obs_H'};
Events = sortrows(Events,'obs_H');

mymap = [
    hex2rgb('1d55c7')
    hex2rgb('1d85f1')
    hex2rgb('5ca9ff')
    hex2rgb('9cd3ff')
    ];

figure('units', 'normalized', 'outerposition', [0 0 1 1]);

scatter(Events.U,Events.obs_H,150,Events.F/1000,'filled');
cb = colorbar;
colormap(mymap)
ylabel(cb,'Fetch [km]','FontSize',16,'Rotation',270)
hold on;
grid on;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THEORY

avg_fetch = mean(Events.F);
std_fetch = std(Events.F);
p_fetch = avg_fetch + std_fetch;
m_fetch = avg_fetch - std_fetch;

% PIERSON-MOSKOWITZ
test_speeds = 1:3:20;
PM = (0.22.*(test_speeds).^2)./9.81;
% JONSWAP
JS = 4.*sqrt((1.67e-7).*((x.^2)./9.81).*avg_fetch);
JS_p = 4.*sqrt((1.67e-7).*((x.^2)./9.81).*p_fetch);
JS_m = 4.*sqrt((1.67e-7).*((x.^2)./9.81).*m_fetch);

plot(x, JS, '--k', 'LineWidth', 3);
plot(x,PM,':k','LineWidth',3)
fill([x, fliplr(x)], [JS_p, fliplr(JS_m)], 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);
fill([x, fliplr(x)], [JS_m, fliplr(JS_p)], 'k', 'FaceAlpha', 0.3, 'EdgeAlpha', 0);


xlabel('$|u|$ [m/s]','FontSize',25,'interpreter','latex')
ylabel('$H_{1/3}$ [m]','FontSize',25,'interpreter','latex')
title('Lake Superior: BUOY 45004')
grid on
box on;
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
ylim([0 max(Events.obs_H)+1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL
planet_to_run = 'Earth';
time_to_run = 1000;   
wind_direction = 0;      
grid_resolution = [20*1000 20*1000];
zDep = 273.5.*ones(12,12);
buoy_loc = [6 6];
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,0,zDep,buoy_loc);
Model.z_data = 3.6;
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);   

% Preallocate cell arrays to store results
avgHsig = cell(1, numel(test_speeds));
pHsig = cell(1, numel(test_speeds));
mHsig = cell(1, numel(test_speeds));
htgrid = cell(1, numel(test_speeds));
E_spec = cell(1, numel(test_speeds));
Cg = cell(1,numel(test_speeds));
wave_height_buoy = NaN(1,numel(test_speeds));

plot(test_speeds, waveh,'--sr','LineWidth',2,'MarkerFaceColor','#c7391d','MarkerSize',15,'MarkerEdgeColor','#c7391d')

for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);

    [avgHsig{i}, ~, E_spec{i}, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
    wave_height_buoy(i) = avgHsig{i}(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT MODEL

    plot(test_speeds(i), avgHsig{i}(end),'--sr','LineWidth',1,'MarkerFaceColor','#c7391d','MarkerSize',15,'MarkerEdgeColor','#c7391d')
    drawnow;

end

% add dashed line between points
plot(test_speeds, wave_height_buoy,'--sr','LineWidth',2,'MarkerFaceColor','#c7391d','MarkerSize',15,'MarkerEdgeColor','#c7391d')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rgb = hex2rgb(hex)
    hex = reshape(hex, [], 6);
    rgb = reshape(sscanf(hex.', '%2x'), [], 3) / 255;
end
