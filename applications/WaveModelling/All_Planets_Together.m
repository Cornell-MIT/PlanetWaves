clc
clear
close all
 
%% PLOT ALL PLANETS TOGETHER FOR SAME DEPTH AND WIND SPEEDS

addpath(genpath(fullfile('..','..','planetwaves')))

save_file = ['AllPlanets_', datestr(datetime("today")), '.mat']; % name of file with saved data

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Parameters
test_speeds = 1:40;
time_to_run = 60*10; % 10 hours  
wind_direction = 0;  

% large, deep basin
grid_resolution = [20*1000 20*1000]; 
zDep = 100.*ones(10,10);
buoy_loc = [5 5];    

% % THRESHOLDS (m/s) (from past tests on this configuration):
% % EARTH                         : 2.2
% % MARS-LOW                      : 1.7
% % MARS-HIGH                     : 1.2
% % TITAN-ONTARIOLACUS            : 0.6
% % TITAN-N2                      : 0.5
% % Kepler-1649b (exo-Venus)      : 5.3
% % LHS-1140b (cold super Earth)  : 2.3
% % 55-Cancri-e (hot super Earth) : 37.1
% 
% % e.g., wavethreshold('Earth') = 2.2

wavethreshold = containers.Map( ...
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
    37.1 ] ...
);

all_planets = {'Earth','Mars-low','Mars-high','Titan-OntarioLacus', 'Titan-N2', 'Kepler-1649-b','LHS-1140-b','55-Cancri-e'};

% wave height vs wind speed per planet
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
    Model.min_freq = 0.01;  % small min freq to allow wave period to grow large for Titan case

    % plot wave height vs time
    figure('Name',['Time Evolution of Waves on ',Planet.name]);
    time_evolve_ax = axes;
    grid on;
    legend('show', 'Location', 'northwest','interpreter','latex');
    title(['Waves on',' ',Planet.name],'interpreter','latex');
    xlabel('model time step [$\Delta$ t]','interpreter','latex')
    ylabel('significant wave height [m]','interpreter','latex')
    hold on;

    % skip non-growth speeds to speed up code
    wind_speeds = test_speeds(test_speeds >= wavethreshold(planet_to_run));
    if ~ismember(wavethreshold(planet_to_run), wind_speeds)
        wind_speeds = [wavethreshold(planet_to_run), wind_speeds];
    end
    u_of_planet{pp} = wind_speeds; % save wind speeds to save file for reference

    time_vs_wave{pp} = cell(numel(wind_speeds),1);    % time evolution of wave height
    wave_height(pp,1:numel(wind_speeds)) = NaN;       % final steady-state height

    for i = 1:numel(wind_speeds)

        Wind.speed = wind_speeds(i);
        Model = calc_cutoff_freq(Planet,Model,Wind);

        % RUN MODEL
        [myHsig{pp,i}, htgrid{pp,i}, ~, ~ , ~ , ~, PeakWaves] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
        
        % peak weighted period
        T_p{pp,i} = PeakWaves.T_weighted;
        % peak wavelength
        L_p{pp,i} = PeakWaves.L;
       
        time_vs_wave{pp}{i} = myHsig{pp,i}; 
        if sum(time_vs_wave{pp}{i}) ~= 0
            wave_height(pp,i) = time_vs_wave{pp}{i}(end);
            plot(time_evolve_ax,1:numel(time_vs_wave{pp}{i}),time_vs_wave{pp}{i},'-','DisplayName',num2str(Wind.speed))
            drawnow
        else
            wave_height(pp,i) = 0;
        end
    end

    save(save_file,"u_of_planet","myHsig","htgrid","T_p")

    % PLOT MODEL
    WAVE_HEIGHT = wave_height(pp,:);
    p1 = plot(sigH_ax,wind_speeds(WAVE_HEIGHT ~= 0), WAVE_HEIGHT(WAVE_HEIGHT ~= 0),'-s','LineWidth',2,'DisplayName',planet_to_run);

    drawnow;

end

legend('show','Location','best')

% write out final results to a table
waveHeight = make_table(all_planets, u_of_planet, wave_height)
% writetable(waveHeight,'WaveHeights.csv','WriteRowNames',true)
