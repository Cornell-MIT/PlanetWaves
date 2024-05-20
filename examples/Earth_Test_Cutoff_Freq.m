clc
clear
close all

planet_to_run = 'Earth';
time_to_run = 10;   
wind_direction = 0;      
grid_resolution = [20*1000 20*1000];
zDep = 273.5.*ones(12,12);
buoy_loc = [6 6];
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,0,zDep,buoy_loc);
Model.z_data = 3.6;
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);   
Wind.speed = 5;
test_cutoff_freq = 10:20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL

figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1])


for j = 1:numel(test_cutoff_freq)
    t = (j-1)/(11-1);
    mycolor = (1-t)*[1,1,0] + t*[0.5,0,0];

    Model.cutoff_freq = test_cutoff_freq(j);

    [myHsig{j}, ~,~, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  

    plot(myHsig{j},'LineWidth',2,'DisplayName',strcat('cutoff freq index =',num2str(test_cutoff_freq(j))),'Color',mycolor)
    hold on
    drawnow;
end

legend('show','Location','eastoutside')
xlabel('model time step [$\Delta$ t]','interpreter','latex')
ylabel('significant wave height [m]','interpreter','latex')
title(strcat('Cut-off Frequency on ',Planet.name))