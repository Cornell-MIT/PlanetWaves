clc
clear
close all

% PLOT ALL PLANETS TOGETHER FOR SAME DEPTH AND WIND SPEEDS

addpath(fullfile('..','..','planetwaves'))  
addpath(fullfile('..','..','planetwaves/pre_analysis/'))  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL
test_speeds = 1:40;
time_to_run = 60*10; % 10 hours  
wind_direction = 0;  
grid_resolution = [20*1000 20*1000];
zDep = 100.*ones(10,10);
buoy_loc = [5 5];    

% THRESHOLDS (m/s) (from past tests on this configuration):
% EARTH                         : 2.2
% MARS-LOW                      : 1.7
% MARS-HIGH                     : 1.2
% TITAN-ONTARIOLACUS            : 0.6
% TITAN-N2                      : 0.5
% Kepler-1649b (exo-Venus)      : 5.3
% LHS-1140b (cold super Earth)  : 2.3
% 55-Cancri-e (hot super Earth) : 36.4

% e.g., wavethreshold('Earth') = 2.2


Planet.wavethreshold = containers.Map( ...
    {'Earth', ...
    'Mars-low', ...
    'Mars-high', ...
    'Titan-OntarioLacus', ...
    'Titan-N2', ...
     'Kepler-1649-b', ...
     'LHS-1140-b', ...
     '55-Cancri-e'}, ...
    [ 2.2, ...
    1.7, ...
    1.2, ...
    0.6, ...
    0.5, ...
    5.3, ...
    2.3, ...
    36.4 ] ...
);

all_planets = {'Earth','Mars-low','Mars-high','Titan-OntarioLacus', 'Titan-N2', 'Kepler-1649-b','LHS-1140-b','55-Cancri-e'};


figure('Name','Sig Wave Heights');
sigH_ax = axes;
xlabel('$|u|$ [m/s]','FontSize',25,'interpreter','latex')
ylabel('$H_{1/3}$ [m]','FontSize',25,'interpreter','latex')
grid on
box on;
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
hold on;


for pp = 1:numel(all_planets)

    planet_to_run = all_planets{pp};

    [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
    Model.gridX = grid_resolution(1);                                              
    Model.gridY = grid_resolution(2);   

    figure('Name',['Time Evolution of Waves on ',Planet.name]);
    time_evolve_ax = axes;
    grid on;
    legend('show', 'Location', 'northwest','interpreter','latex');
    title(['Waves on',' ',Planet.name],'interpreter','latex');
    xlabel('model time step [$\Delta$ t]','interpreter','latex')
    ylabel('significant wave height [m]','interpreter','latex')
    hold on;
    

    if strcmp(planet_to_run,'55-Cancrie') % skip non-growth values to run faster
        test_speeds = [35:40];
    end
    
 
    time_vs_wave = NaN(numel(test_speeds),time_to_run);
    for i = 1:numel(test_speeds)
    
        Wind.speed = test_speeds(i);
        Model = calc_cutoff_freq(Planet,Model,Wind);

        [avgHsig, ~, ~, ~, ~, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
        time_vs_wave(i,:) = avgHsig;
        save_avgHsig = avgHsig;
        save_avgHsig(avgHsig==0) = [];
        if sum(avgHsig) ~= 0
            wave_height(pp,i) = save_avgHsig(end);
            plot(time_evolve_ax,1:numel(avgHsig),avgHsig,'-','DisplayName',num2str(Wind.speed))
            drawnow;
            hold on;
            yline(time_evolve_ax,wave_height(pp,i),'--k',num2str(Wind.speed),'DisplayName',['H at ' num2str(Wind.speed)])
            drawnow;
        else
            wave_height(pp,i) = 0;
        end
    
    
    end
    
    % PLOT MODEL
    WAVE_HEIGHT = wave_height(pp,:);
    p1 = plot(sigH_ax,test_speeds(WAVE_HEIGHT ~= 0), WAVE_HEIGHT(WAVE_HEIGHT ~= 0),'-s','LineWidth',2,'DisplayName',planet_to_run);
  
    drawnow;

end

legend('show','Location','best')