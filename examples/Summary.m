clc
clear
close all

% PLOT ALL PLANETS TOGETHER FOR SAME DEPTH AND WIND SPEEDS

addpath(fullfile('..','planetwaves'))  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL
test_speeds = 1:10;
time_to_run = 60*10;  
wind_direction = 0;  
buoy_loc = [6 6];    
grid_resolution = [7.5*1000 7.5*1000];
zDep = 10.*ones(12,12);

all_planets = {'Earth','Mars','Mars-high','Titan-OntarioLacus','Exo-Venus','Titan-N2','LHS-1140b'};
%all_planets = {'Titan'};

figure;
sigH_ax = axes;
xlabel('$|u|$ [m/s]','FontSize',25,'interpreter','latex')
ylabel('$H_{1/3}$ [m]','FontSize',25,'interpreter','latex')
grid on
box on;
xlim([0 max(test_speeds)+1])
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
hold on;


for pp = 1:numel(all_planets)

    planet_to_run = all_planets{pp};

    [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
    Model.gridX = grid_resolution(1);                                              
    Model.gridY = grid_resolution(2);   

    figure;
    time_evolve_ax = axes;
    grid on;
    legend('show', 'Location', 'northwest','interpreter','latex');
    title(['Waves on',' ',Planet.name],'interpreter','latex');
    xlabel('model time step [$\Delta$ t]','interpreter','latex')
    ylabel('significant wave height [m]','interpreter','latex')
    hold on;
    
    time_vs_wave = NaN(numel(test_speeds),time_to_run);
    for i = 1:numel(test_speeds)
    
        Wind.speed = test_speeds(i);
        Model = calc_cutoff_freq(Planet,Model,Wind);
    
        [avgHsig, ~, ~, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
        time_vs_wave(i,:) = avgHsig;
        save_avgHsig = avgHsig;
        save_avgHsig(avgHsig==0) = [];
        if sum(avgHsig) ~= 0
            spectrogram.wave_height(i) = save_avgHsig(end);
            plot(time_evolve_ax,1:numel(avgHsig),avgHsig,'-','DisplayName',num2str(Wind.speed))
            drawnow;
            hold on;
            yline(time_evolve_ax,spectrogram.wave_height(i),'--k',num2str(Wind.speed),'DisplayName',['H at ' num2str(Wind.speed)])
            drawnow;
        else
            spectrogram.wave_height(i) = 0;
        end
        spectrogram.wind(i) = test_speeds(i);
    
    
    end
    
    % PLOT MODEL
    p1 = plot(sigH_ax,test_speeds(spectrogram.wave_height ~=0), spectrogram.wave_height(spectrogram.wave_height~=0),'-s','LineWidth',2,'DisplayName',planet_to_run);
    drawnow;

end

legend('show','Location','best')