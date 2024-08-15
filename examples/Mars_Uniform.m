clc
clear
close all

% COMPARE THE EFFECT OF CHANGING THE ATMOSPHERIC PRESSURE ON LAKES AT MARS

addpath(fullfile('..','planetwaves'))  

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL
test_speeds = [1.2 1.6 2:10];
planet_to_run = 'Mars';
time_to_run = 60*10;   
wind_direction = 0;      
grid_resolution = [7.5*1000 7.5*1000];
zDep = 10.*ones(12,12);
buoy_loc = [6 6];
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
Model.z_data = 10;
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);   


figure
surf(zDep);
view(2)
new_xtick = get(gca, 'XTick')*(Model.gridX)/1000;
new_ytick = get(gca, 'YTick')*(Model.gridY)/1000;
colorbar
set(gca, 'XTick',  get(gca, 'XTick'), 'XTickLabel', arrayfun(@(x) sprintf('%.2f', x), new_xtick, 'UniformOutput', false));
set(gca, 'YTick',  get(gca, 'YTick'), 'YTickLabel', arrayfun(@(y) sprintf('%.2f', y), new_ytick, 'UniformOutput', false));
xlabel('longitude [km]')
ylabel('latitude [km]')
title('model input')

figure;
sigH_ax = axes;
xlabel('$|u|$ [m/s]','FontSize',25,'interpreter','latex')
ylabel('$H_{1/3}$ [m]','FontSize',25,'interpreter','latex')
title('Mars: Jezero Crater Lake')
grid on
box on;
xlim([0 max(test_speeds)+1])
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')
hold on;

figure;
time_evolve_ax = axes;
grid on;
legend('show', 'Location', 'northwest','interpreter','latex');
title(['Waves on',' ',Planet.name,' at ',num2str(Planet.surface_press/1000), ' kPa'],'interpreter','latex');
xlabel('model time step [$\Delta$ t]','interpreter','latex')
ylabel('significant wave height [m]','interpreter','latex')
hold on;

time_vs_wave = NaN(numel(test_speeds),time_to_run);

for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);
    Model = calc_cutoff_freq(Planet,Model,Wind);

    [avgHsig, ~, E_spec, ~, ~,~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
    time_vs_wave(i,:) = avgHsig;
    save_avgHsig = avgHsig;
    save_avgHsig(avgHsig==0) = [];
    save_E_spec = E_spec(~cellfun('isempty',E_spec));
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
    spectrogram.energy{i} = save_E_spec{end};
    spectrogram.wind(i) = test_speeds(i);

    save('LongMarsRunLowPres.mat','spectrogram')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT MODEL
    if spectrogram.wave_height(i) > 0
        plot(sigH_ax,test_speeds(i), spectrogram.wave_height(i),'--s','LineWidth',1,'MarkerFaceColor','#FF6F59','MarkerSize',15,'MarkerEdgeColor','#FF6F59','Color','#FF6F59')
        drawnow;
    end

end



% add dashed line between points
p1 = plot(sigH_ax,test_speeds, spectrogram.wave_height,'--s','LineWidth',2,'MarkerFaceColor','#FF6F59','MarkerSize',15,'MarkerEdgeColor','#FF6F59','MarkerEdgeColor','#FF6F59','Color','#FF6F59');



Planet.surface_press = 4*Planet.surface_press;

figure;
time_evolve_ax2 = axes;
grid on;
legend('show', 'Location', 'northwest','interpreter','latex');
title(['Waves on',' ',Planet.name,' at ',num2str(Planet.surface_press/1000), ' kPa'],'interpreter','latex');
xlabel('model time step [$\Delta$ t]','interpreter','latex')
ylabel('significant wave height [m]','interpreter','latex')
hold on;


time_vs_wave2 = NaN(numel(test_speeds),time_to_run);

for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);
    Model = calc_cutoff_freq(Planet,Model,Wind);

    [avgHsig, ~, E_spec, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc); 
    time_vs_wave2(i,:) = avgHsig;
    save_avgHsig = avgHsig;
    save_avgHsig(avgHsig==0) = [];
    save_E_spec = E_spec(~cellfun('isempty',E_spec));
    if sum(avgHsig) ~= 0
        spectrogram2.wave_height(i) = save_avgHsig(end);
        plot(time_evolve_ax2,1:numel(avgHsig),avgHsig,'-','DisplayName',num2str(Wind.speed))
        drawnow;
        hold on;
        yline(time_evolve_ax2,spectrogram2.wave_height(i),'--k',num2str(Wind.speed),'DisplayName',['H at ' num2str(Wind.speed)])
        drawnow;
    else
        spectrogram2.wave_height(i) = 0;
    end
    spectrogram2.energy{i} = save_E_spec{end};
    spectrogram2.wind(i) = test_speeds(i);
    save('LongMarsRunLowPres.mat','spectrogram2')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT MODEL
    if spectrogram2.wave_height(i) > 0
        plot(sigH_ax,test_speeds(i), spectrogram2.wave_height(i),'--s','LineWidth',1,'MarkerFaceColor','#2B3A67','MarkerSize',15,'MarkerEdgeColor','#2B3A67','Color','#2B3A67')
        drawnow;
    end

end


% add dashed line between points
p2 = plot(sigH_ax,test_speeds, spectrogram2.wave_height,'--s','LineWidth',2,'MarkerFaceColor','#2B3A67','MarkerSize',15,'MarkerEdgeColor','#2B3A67','MarkerEdgeColor','#2B3A67','Color','#2B3A67');

legend([p1 p2],'50 kPA','200 kPa','Location','best')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rgb = hex2rgb(hex)
    hex = reshape(hex, [], 6);
    rgb = reshape(sscanf(hex.', '%2x'), [], 3) / 255;
end