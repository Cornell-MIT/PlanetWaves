clc
clear
close all

addpath(fullfile('..','planetwaves'))  

planet_to_run = 'Exo-CO2';
test_speeds = [0.5 1.0 2.0 5.0 8.0 11.0 14.0 17.0 20.0]; % initiates at around 0.8 m/s
time_to_run = 500;   
wind_direction = 0;      
grid_resolution = [20*1000 20*1000];
zDep = 273.5.*ones(12,12);
buoy_loc = [6 6];

[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);

avgHsig = cell(1, numel(test_speeds));
E_spec = cell(1, numel(test_speeds));
wave_height_buoy = NaN(1,numel(test_speeds));
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);                                                                                                 


for i = 1:numel(test_speeds)

    t = (i-1)/(numel(test_speeds)-1);
    mycolor = (1-t)*[0,1,1] + t*[0,0,1];

    Wind.speed = test_speeds(i);

    [avgHsig{i}, ~, E_spec{i}, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
    wave_height_buoy(i) = avgHsig{i}(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT MODEL
    if avgHsig{i}(end) > 0
        plot(test_speeds(i), avgHsig{i}(end),'sb','MarkerFaceColor','b','MarkerSize',15,'MarkerEdgeColor','b')
    end
    hold on;
    drawnow;

end
grid on;
%legend('show','Location','eastoutside')
xlabel('|u| [m/s]','interpreter','latex')
ylabel('significant wave height [m]','interpreter','latex')
% add dashed line between points
plot(test_speeds, wave_height_buoy,'--b','LineWidth',2)

planet_to_run = 'Earth';
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);

avgHsig2 = cell(1, numel(test_speeds));
E_spec2 = cell(1, numel(test_speeds));
wave_height_buoy = NaN(1,numel(test_speeds));
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);                                                                                                 


for i = 1:numel(test_speeds)

    t = (i-1)/(numel(test_speeds)-1);
    mycolor = (1-t)*[0,1,1] + t*[0,0,1];

    Wind.speed = test_speeds(i);

    [avgHsig2{i}, ~, E_spec2{i}, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
    wave_height_buoy(i) = avgHsig2{i}(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT MODEL
    if avgHsig2{i}(end) > 0
        plot(test_speeds(i), avgHsig2{i}(end),'sg','MarkerFaceColor','g','MarkerSize',15,'MarkerEdgeColor','g')
    end
    drawnow;

end

plot(test_speeds, wave_height_buoy,'--sg','LineWidth',2)

planet_to_run = 'Titan';
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);

avgHsig3 = cell(1, numel(test_speeds));
E_spec3 = cell(1, numel(test_speeds));
wave_height_buoy = NaN(1,numel(test_speeds));
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);                                                                                                 


for i = 1:numel(test_speeds)

    t = (i-1)/(numel(test_speeds)-1);
    mycolor = (1-t)*[0,1,1] + t*[0,0,1];

    Wind.speed = test_speeds(i);

    [avgHsig3{i}, ~, E_spec3{i}, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
    wave_height_buoy(i) = avgHsig3{i}(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT MODEL
    if avgHsig3{i}(end) > 0
        plot(test_speeds(i), avgHsig3{i}(end),'sr','MarkerFaceColor','r','MarkerSize',15,'MarkerEdgeColor','r')
    end
    drawnow;

end

plot(test_speeds, wave_height_buoy,'--sr','LineWidth',2)
