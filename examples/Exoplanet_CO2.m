clc
clear
close all

addpath(fullfile('..','planetwaves'))  

planet_to_run = 'Exo-CO2';
test_speeds = 0.1:0.1:10; % initiates at around 0.8 m/s
time_to_run = 10;   
wind_direction = 0;      
grid_resolution = [20*1000 20*1000];
zDep = 273.5.*ones(12,12);
buoy_loc = [6 6];

[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,0,zDep,buoy_loc);

avgHsig = cell(1, numel(test_speeds));
E_spec = cell(1, numel(test_speeds));
wave_height_buoy = NaN(1,numel(test_speeds));
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);   
Model.cutoff_freq = 25;
Model.mindelt = 0.1;
Model.maxdelt = 500.0; 
for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);

    [avgHsig{i}, ~, E_spec{i}, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
    wave_height_buoy(i) = avgHsig{i}(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT MODEL
    plot(avgHsig{i})
    hold on;
    %plot(test_speeds(i), avgHsig{i}(end),'--sr','LineWidth',1,'MarkerFaceColor','#c7391d','MarkerSize',15,'MarkerEdgeColor','#c7391d')
    drawnow;

end

% add dashed line between points
%plot(test_speeds, wave_height_buoy,'--sr','LineWidth',2,'MarkerFaceColor','#c7391d','MarkerSize',15,'MarkerEdgeColor','#c7391d')