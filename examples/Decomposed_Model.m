clc
clear
close all

% SHOW DECOMPOSITION OF MODEL 

addpath(fullfile('..','planetwaves'))  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL
test_speeds = 3;
planet_to_run = 'Earth';
time_to_run = 60*10;   
wind_direction = 0;      
grid_resolution = [20*1000 20*1000];
zDep = 273.5.*ones(12,12);
buoy_loc = [6 6];
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,0,zDep,buoy_loc);
Model.z_data = 3.6;
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);
Model.tolH = NaN;
Etc.showplots = 1;

for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);

    [avgHsig, ~, ~, ~, ~,~,~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 

 
end
print('-vectpr', '-dpdf', 'Decomposed_Model')
x = 1:numel(avgHsig);
y = avgHsig;
z = zeros(size(avgHsig));
col = flip(autumn(time_to_run),1); % yellow -> red, with 61 colors (for 61 lines)
figure;
for ii = 1:20:numel(avgHsig)
    plot(x(ii),y(ii),'.','Color',col(ii,:),'MarkerSize',30)
    hold on;
end
grid on;
xlabel('model time step [$\Delta$ t]','interpreter','latex')
ylabel('significant wave height [m]','interpreter','latex')
% print('-vector', '-dpdf', 'Decomposed_Waveheight')
