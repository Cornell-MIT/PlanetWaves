clc
clear
close all

addpath(fullfile('..','planetwaves'))  
addpath(fullfile('..','planetwaves','pre_analysis'))  

% MAKE PIERSON-MOSKOWITZ-ish CURVE FOR CHANGING PLANETARY CONDITIONS FOR HSIG AND UMIN

planet_to_run = 'Earth';
test_speeds = 10;
time_to_run = 60*10;  
wind_direction = pi/2;  
buoy_loc = [3 40];    
grid_resolution = [0.5*1000 1*1000];
zDep = 100.*ones(80,6);

[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);   


% (1)
%Earth = Planet.kgmolwt;
rat = [0.25 0.5 1 2 4];
%var = Earth.*rat;
%test_speeds = test_speeds.*rat;


figure;
time_evolve_ax = axes;
grid on;
legend('show', 'Location', 'northwest','interpreter','latex');
title(['Waves on',' ',Planet.name,' at ',num2str(Planet.surface_press/1000), ' kPa'],'interpreter','latex');
xlabel('model time step [$\Delta$ t]','interpreter','latex')
ylabel('significant wave height [m]','interpreter','latex')
hold on;


for i = 1:numel(test_speeds)

                Wind.speed = test_speeds(i);
% (2)
                %Planet.kgmolwt = var(i);
                Model = calc_cutoff_freq(Planet,Model,Wind);

                make_input_map(Planet,Model,Wind)
                [myHsig,htgrid, ~, ~ , ~ , ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
                plot(time_evolve_ax,1:numel(myHsig),myHsig,'-','DisplayName',num2str(Wind.speed))
                Hsig(i) = myHsig(end);

end

% my_rat = Hsig./Hsig(rat == 1)
% 
% % fit power law
% X = rat; 
% Y = my_rat;
logX = log(X);
logY = log(Y);
fitResult = polyfit(logX, logY, 1);
n = fitResult(1);
fprintf('The best fit exponent n is: %f\n', n);





