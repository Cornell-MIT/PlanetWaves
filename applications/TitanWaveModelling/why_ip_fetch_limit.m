clc
clear
close all

% Why is there an inflection point in the fetch limit for increasing wind speeds?

addpath(fullfile('..','..','planetwaves'))  
addpath(fullfile('..','..','planetwaves/pre_analysis/'))  


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL
test_speeds = 40;
time_to_run = 60*10;  
wind_direction = 0;  
buoy_loc = [5 5];    
grid_resolution = [10*1000 10*1000];
zDep = 100.*ones(10,9);

all_planets = {'Titan-OntarioLacus'};




for pp = 1:numel(all_planets)

    planet_to_run = all_planets{pp};

    [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
    Model.gridX = grid_resolution(1);                                              
    Model.gridY = grid_resolution(2);   
    Etc.showplots = 1;

 
    time_vs_wave = NaN(numel(test_speeds),time_to_run);
    for i = 1:numel(test_speeds)
    
        Wind.speed = test_speeds(i);
        Model = calc_cutoff_freq(Planet,Model,Wind);
        disp(Model.cutoff_freq(i))

        [avgHsig, ~, ~, ~, ~, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
        time_vs_wave(i,:) = avgHsig;
        save_avgHsig = avgHsig;
        save_avgHsig(avgHsig==0) = [];
        if sum(avgHsig) ~= 0
            wave_height(pp,i) = save_avgHsig(end);

        else
            wave_height(pp,i) = 0;
        end
    
    
    end
    


end

