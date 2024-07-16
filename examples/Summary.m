clc
clear
close all

% plot all planets together for same depth profile

addpath(fullfile('..','planetwaves'))  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL
test_speeds = 0:10;
time_to_run = 50;   
    
grid_resolution = [7.5*1000 7.5*1000];
zDep = 10.*ones(12,12);

figure
surf(zDep);
view(2)
new_xtick = get(gca, 'XTick')*(grid_resolution(1))/1000;
new_ytick = get(gca, 'YTick')*(grid_resolution(2))/1000;
colorbar
set(gca, 'XTick',  get(gca, 'XTick'), 'XTickLabel', arrayfun(@(x) sprintf('%.2f', x), new_xtick, 'UniformOutput', false));
set(gca, 'YTick',  get(gca, 'YTick'), 'YTickLabel', arrayfun(@(y) sprintf('%.2f', y), new_ytick, 'UniformOutput', false));
xlabel('longitude [km]')
ylabel('latitude [km]')
title('model input bathymetry')

%all_planets = {'Earth','Mars','Titan','Exo-Venus'};
all_planets = {'Titan'};
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

    wind_direction = 0;  
    buoy_loc = [6 6];
    [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
    Model.z_data = 10;
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
    p1 = plot(sigH_ax,test_speeds, spectrogram.wave_height,'--s','LineWidth',2,'DisplayName',planet_to_run);
    drawnow


end

legend('show','Location','best')


