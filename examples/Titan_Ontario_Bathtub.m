clc
clear
close all

% ONTARIO LACUS WITH BATHTUB BATHYMETRY

addpath(fullfile('..','planetwaves'))  
addpath(fullfile('..','planetwaves','pre_analysis'))
addpath(fullfile('..','data','Titan','TitanLakes','Bathymetries','bathtub_bathy'))
load('..\data\Titan\TitanLakes\Bathymetries\SAR_bathy_cleaned\ol_main_basin.mat','smoothed_ol');

zDep = smoothed_ol;
zDep = imrotate(zDep,180);
zDep = imrotate(zDep,-90);
zDep(:,80:end) = [];
zDep(1:95,:) = [];
zDep(90:end,:) = [];


% MODEL INPUTS
planet_to_run = 'Titan';
buoy_loc = [60 30];                                                        % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
test_speeds = 1;                                                     % wind speed
time_to_run = 720;                                                          % time to run model
wind_direction = pi/2;                                                        % wind direction


buoy_loc = [400 800];                                                      % grid location [x,y]
grid_resolution = [180 180];                                       % pixel width and pixel height [m]
test_speeds = 1:3;                                                 % wind speed
time_to_run = 100;                                                        % time to run model
wind_direction = 0;                                                        % wind direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL MODEL
% degrade depth profile so model doesnt take as long to run
[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,1);
% populate model classes
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
% update grid resolution
Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);                                               
          
figure;
imagesc(zDep)
hold on;
plot(buoy_loc(1),buoy_loc(2),'or','MarkerFaceColor','r')
colorbar;
hold off;
title('input bathymetry')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL

% Preallocate cell arrays to store results
myHsig = cell(1, numel(test_speeds));
htgrid = cell(1, numel(test_speeds));
E_spec = cell(1, numel(test_speeds));
Cg = cell(1,numel(test_speeds));

for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);
    
    [myHsig{i}, htgrid{i}, ~, ~, ~,wave_age{i}] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PLOT RESULTS

figure
for i = 1:numel(test_speeds)
    plot(myHsig{i}, '-', 'LineWidth', 3, 'DisplayName', num2str(test_speeds(i)))
    hold on
end
grid on;
legend('show', 'Location', 'northwest','interpreter','latex');
title(['Waves on',' ',Planet.name],'interpreter','latex');
xlabel('model time step [$\Delta$ t]','interpreter','latex')
ylabel('significant wave height [m]','interpreter','latex')


for k = 1:numel(test_speeds)
    buoy_waves(k) = htgrid{k}{end}(Model.long,Model.lat);
end

figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1])
for speed = 1:numel(test_speeds)
    
    subplot(1,3,3)
    plot(test_speeds,buoy_waves,'-sb','LineWidth',1,'MarkerFaceColor','b')
    hold on
    plot(test_speeds(speed),htgrid{speed}{end}(Model.long,Model.lat),'-sr','LineWidth',1,'MarkerFaceColor','r')
    hold off;
    xlabel('u [m/s]','interpreter','latex')
    ylabel('H_{sig} [m]','interpreter','latex')
    grid on;
    
    plot_grid = htgrid{speed}{end}';
    plot_grid(isnan(plot_grid)) = 0;
    plot_alpha_data = ones(size(plot_grid));
    plot_alpha_data(plot_grid==0) = 0;

    ax1 = subplot(1,3,[1,2]);
    h1 = imagesc(plot_grid);
    colormap cool
    xlabel('longitude [km]')
    ylabel('latitude [km]')
    title(sprintf('u = %i m/s',test_speeds(speed)))
    c1 = colorbar;
    c1.Label.String = 'H_{sig} [m]';
    clim([0 2])
    set(h1, 'AlphaData', plot_alpha_data);
    hold on;
    
    contour(Model.bathy_map,min(min(Model.bathy_map)):10:max(max(Model.bathy_map)),'-k','LineWidth',2)

    grid on
    new_xtick = get(gca, 'XTick')*(Model.gridX)/1000;
    new_ytick = get(gca, 'YTick')*(Model.gridY)/1000;
    set(gca, 'XTick',  get(gca, 'XTick'), 'XTickLabel', arrayfun(@(x) sprintf('%d', x), new_xtick, 'UniformOutput', false));
    set(gca, 'YTick',  get(gca, 'YTick'), 'YTickLabel', arrayfun(@(y) sprintf('%d', y), new_ytick, 'UniformOutput', false));
    
    
    [wx,wy] = pol2cart(Wind.dir,1);
    plot(Model.long,Model.lat,'pentagram','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20)

    for i = 1:Model.LonDim
        for j = 1:Model.LatDim
            if Model.bathy_map(j,i) <=0
                quiver(i, j, (speed/5)*wx, (speed/5)*wy, 'r', 'MaxHeadSize', 1);
            end
        end
    end
    set(ax1,'Ydir','reverse')
    
    drawnow
    if speed == 1
        gif('Ontario_Titan_Bathtub.gif','DelayTime',1,'overwrite',true)
    else
        gif
    end
end



figure('units','normalized','outerposition',[0 0 1 1])
plot(test_speeds,buoy_waves,'-sb','LineWidth',1,'MarkerFaceColor','b')
%ylim([0 3])
xlabel('u [m/s]','Interpreter','latex')
ylabel('H_{sig} [m]','Interpreter','tex')
grid on;


maturity = wave_age{1};
maturity(maturity>=0.83) = 1;
maturity(maturity<0.83) = 0;

figure
h2 = imagesc(maturity);
cRange = caxis; 
hold on;
contour(Model.bathy_map,min(min(Model.bathy_map)):10:max(max(Model.bathy_map)),'-k','LineWidth',2)
caxis(cRange);
title('Maturity of Waves')
set(h2, 'AlphaData', plot_alpha_data);
