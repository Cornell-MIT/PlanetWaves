clc
clear
close all

addpath(fullfile('..','planetwaves'))                                                     % location of wave model files

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETERS:
% (1) PLANET CONDITIONS
%   (a) TITAN
Titan.rho_liquid = 540;                                                    % Hydrocarbon Liquid density [kg/m3]
Titan.nu_liquid = 3e-7;                                                    % Hydrocarbon Liquid Viscosity [m2/s]
Titan.nua = 0.0126/1e4;                                                    % Titan atmospheric gas viscosity [m2/s]
Titan.gravity = 1.352;                                                     % Titan Gravity [m/s2]
Titan.surface_temp = 92;                                                   % Titan Surface Temperature [K]
Titan.surface_press = 1.5*101300;                                          % Titan Surface Pressure [Pa]
Titan.surface_tension = 0.018;                                             % Hydrocarbon Liquid Surface Tension [N/m]
Titan.name = 'Titan';
%   (b) EARTH
Earth.rho_liquid = 997;                                                    % Water Liqid Density [kg/m3]       
Earth.nu_liquid = 1.143e-6;                                                % Water Liquid Viscosity [m2/s]
Earth.nua = 1.48e-5;                                                       % Earth atmospheric gas viscosity [m2/s]
Earth.gravity = 9.81;                                                      % Earth Gravity [m/s2]
Earth.surface_temp = 288;                                                  % Earth Surface Temperature [K]
Earth.surface_press = 1*101300;                                            % Earth Surface Pressure [Pa]
Earth.surface_tension = 0.072;                                             % Water Liquid Surface Tension [N/m]
Earth.name = 'Earth';
% (2a) MODEL GEOMETRY
Model.LonDim = 10;                                                         % Number of Grid Cells in X-Dimension (col count)
Model.LatDim = 10;                                                         % Number of Grid Cells in Y-Dimension (row count)
Model.Fdim = 35;                                                           % Number of Frequency bins
Model.Dirdim = 72;                                                         % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.long = 5;                                                            % longitude grid point for sampling during plotting
Model.lat = 5;                                                             % latitude grid point for sampling during plotting
Model.gridX = 40*1000;                                                     % Grid size in X-dimension [m]
Model.gridY = 40*1000;                                                     % Grid size in Y-dimension [m]
Model.mindelt = 0.0001;                                                    % minimum time step
Model.maxdelt = 2000.0;                                                    % maximum time step
Model.time_step = 100;                                                     % Maximum Size of time step [s] -- if set too low can lead to numerical ringing
Model.num_time_steps = 200;                                                % Length of model run (in terms of # of time steps)
Model.tolH = NaN;                                                          % tolerance threshold for maturity
Model.cutoff_freq = 15;                                                    % cutoff frequency bin from diagnostic to advection -- if set too low can lead to numerical ringing
Model.min_freq = 0.05;                                                     % minimum frequency to model
Model.max_freq = 35;                                                       % maximum frequency to model
% (2b) BUOY SPECIFC:
% STATION 45012:
% https://www.ndbc.noaa.gov/station_page.php?station=45012
Model.z_data = 3.6;                                                        % elevation of wind measurement [m]
Model.bathy_map = 273.6.*ones(Model.LonDim,Model.LatDim);                  % depth of water column beneath buoy [m]
% (2c) TUNING PARAMETERS
Model.tune_A1 = 0.11;                                                      % wind sea (eq. 12 Donelan+2012)
Model.tune_mss_fac = 360;
Model.tune_Sdt_fac = 0.001;
Model.tune_Sbf_fac = 0.002;
Model.tune_cotharg = 0.2;
Model.tune_n = 2.4;
% (3) NEAR-SURFACE WIND CONDITIONS
test_speeds = [3 8 10];                                                    % magnitude of incoming wind [m/s] [e.g.
Wind.dir = (3*pi)/2;                                                       % direction of incoming wind [radians] where 0 deg corresponds to winds traveling from left to right and increases clockwise
% (4) Unidirectional currents
Uniflow.East = 0;                                                          % eastward unidirectional current [m/s]
Uniflow.North = 0;                                                         % northward unidirectional current [m/s]
% (5) HOUSEKEEPING
Etc.showplots = 0;
Etc.savedata = 1;

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% RUN THE MODEL
planet_to_run = Earth;
myHsig = NaN(numel(test_speeds),Model.num_time_steps);
for i = 1:numel(test_speeds)
	Wind.speed  = test_speeds(i);
	[myHsig(i,:),htgrid{i},~,~] = makeWaves(planet_to_run,Model,Wind,Uniflow,Etc);
  if i == 1
      figure('units','normalized','outerposition',[0 0 1 1])
  end
  plot(myHsig(i,:),'-','LineWidth',3,'DisplayName', num2str(Wind.speed))
  hold on
  drawnow
end
grid on;
legend('show', 'Location', 'northwest','interpreter','latex');
title(['Waves on',' ',planet_to_run.name],'interpreter','latex');
xlabel('model time step [$\Delta$ t]','interpreter','latex')
ylabel('significant wave height [m]','interpreter','latex')
disp('run finished')

% plot the results
for k = 1:numel(test_speeds)
    buoy_waves(k) = htgrid{k}{end}(Model.long,Model.lat);
end

figure('units','normalized','outerposition',[0 0 1 1])
for speed = 1:numel(test_speeds)
    
    subplot(1,3,3)
    plot(test_speeds,buoy_waves,'-sb','LineWidth',1,'MarkerFaceColor','b')
    hold on
    plot(test_speeds(speed),htgrid{speed}{end}(Model.long,Model.lat),'-sr','LineWidth',1,'MarkerFaceColor','r')
    hold off;
    xlabel('u [m/s]')
    ylabel('H_{sig} [m]')
    grid on;
    
    plot_grid = htgrid{speed}{end}';
    plot_grid(isnan(plot_grid)) = 0;
    plot_alpha_data = ones(size(plot_grid));
    plot_alpha_data(plot_grid==0) = 0;

    ax1 = subplot(1,3,[1,2]);
    h1 = imagesc(plot_grid);
    colormap linspecer
    xlabel('longitude [km]')
    ylabel('latitude [km]')
    title(sprintf('u = %i m/s',test_speeds(speed)))
    c1 = colorbar;
    c1.Label.String = 'H_{sig} [m]';
    set(h1, 'AlphaData', plot_alpha_data);
    hold on;
    
    contour(plot_grid,'-k','LineWidth',2)

    grid on
    new_xtick = get(gca, 'XTick')*(Model.gridX)/1000;
    new_ytick = get(gca, 'YTick')*(Model.gridY)/1000;
    set(gca, 'XTick',  get(gca, 'XTick'), 'XTickLabel', arrayfun(@(x) sprintf('%d', x), new_xtick, 'UniformOutput', false));
    set(gca, 'YTick',  get(gca, 'YTick'), 'YTickLabel', arrayfun(@(y) sprintf('%d', y), new_ytick, 'UniformOutput', false));
    
    
    [wx,wy] = pol2cart(Wind.dir,1);
    plot(Model.long,Model.lat,'pentagram','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20)

    for i = 1:Model.LonDim
        for j = 1:Model.LatDim
            quiver(i, j, (speed/5)*wx, (speed/5)*wy, 'k', 'MaxHeadSize', 1);
        end
    end
    set(ax1,'Ydir','reverse')
    drawnow
    if speed == 1
        gif('example.gif','DelayTime',1)
    else
        gif
    end
end


