clc
clear
close all
% from find_fetch.py
load('.\EarthAnalysis\GreatLakes\LakeData\LakeSuperior_cleaned.mat')
LS = squeeze(LS);
LS_orig = -LS;
resizeFactor = 0.002;
% from find_fetch.py
blon = 1729;
blat = 6618;
gridcellsizeX = 4542.948547909539*(1/resizeFactor);
gridcellsizeY = 92.66280063299297*(1/resizeFactor);
pos =  [blat, blon];
pos_orig = pos;
pos = round(pos * resizeFactor);
%LS(LS == 0) = NaN;
LS = imresize(LS, resizeFactor, "bilinear");
alphaData = ones(size(LS));
alphaData(LS==0) = 0;
LS = -LS;
%figure; imagesc(-LS,'AlphaData', alphaData); hold on; plot(position(1),position(2),'ok','MarkerFaceColor','k'); colorbar; colormap cool;
% figure;
% subplot(2,1,1)
% surf(LS_orig,'edgecolor','none'); view(2); hold on; plot3(pos_orig(1),pos_orig(2),1e4,'or','MarkerFaceColor','r'); colormap cool; colorbar; title('original resolution')
% subplot(2,1,2)
% surf(LS,'edgecolor','none'); view(2); hold on; plot3(pos(1),pos(2),1e4,'or','MarkerFaceColor','r'); colormap cool; colorbar; title('degraded resolution')

size_lake = size(LS);
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
Model.m = 10;%size_lake(2);                                                    % Number of Grid Cells in X-Dimension
Model.n = 10;%size_lake(1);                                                    % Number of Grid Cells in Y-Dimension
Model.o = 35;                                                              % Number of Frequency bins
Model.p = 72;                                                              % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.long = pos(1);                                                       % longitude grid point for sampling during plotting
Model.lat = pos(2);                                                        % latitude grid point for sampling during plotting
Model.gridX = gridcellsizeX;                                                     % Grid size in X-dimension [m]
Model.gridY = gridcellsizeY;                                                     % Grid size in Y-dimension [m]
Model.mindelt = 0.0001;                                                    % minimum time step
Model.maxdelt = 2000.0;                                                    % maximum time step
Model.time_step = 100;                                                     % Maximum Size of time step [s] -- if set too low can lead to numerical ringing
Model.num_time_steps = 1000;                                                % Length of model run (in terms of # of time steps)
Model.tolH = NaN;                                                          % tolerance threshold for maturity
Model.cutoff_freq = 15;                                                    % cutoff frequency bin from diagnostic to advection -- if set too low can lead to numerical ringing
Model.min_freq = 0.05;                                                     % minimum frequency to model
Model.max_freq = 35;                                                       % maximum frequency to model
% (2b) BUOY SPECIFC:
% STATION 45012:
% https://www.ndbc.noaa.gov/station_page.php?station=45012
Model.z_data = 3.6;                                                         % elevation of wind measurement [m]
Model.bathy_map = 273.6.*ones(Model.m,Model.n);                             % depth of water column beneath buoy [m]
%Model.bathy_map = LS;                                                      % bathymetry of model basin [m]
% (2c) TUNING PARAMETERS
Model.tune_A1 = 0.11;                                                       % wind sea (eq. 12 Donelan+2012)
Model.tune_mss_fac = 360;
Model.tune_Sdt_fac = 0.001;
Model.tune_Sbf_fac = 0.002;
Model.tune_cotharg = 0.2;
Model.tune_n = 2.4;
% (3) NEAR-SURFACE WIND CONDITIONS
test_speeds = 15;%1:2:20;                                                      % magnitude of incoming wind [m/s] [e.g.
Wind.dir = 0;                                                              % direction of incoming wind [radians]
% (4) Unidirectional currents
Uniflow.East = 0;                                                          % eastward unidirectional current [m/s]
Uniflow.North = 0;                                                         % northward unidirectional current [m/s]
% (5) HOUSEKEEPING
Etc.showplots = 0;
Etc.savedata = 0;
Etc.showlog = 0;
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% RUN THE MODEL
planet_to_run = Earth;
myHsig = NaN(numel(test_speeds),Model.num_time_steps);
for i = 1:numel(test_speeds)
	Wind.speed  = test_speeds(i);
	[myHsig(i,:),htgrid{i},~,~] = makeWaves(planet_to_run,Model,Wind,Uniflow,Etc);
  ht_sig(i) = htgrid{k}{end}(Model.long,Model.lat);
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
figure;
plot(test_speeds(i),ht_sig(i),'--sr','LineWidth',1,'MarkerFaceColor','r')
xlabel('$|u|$ [m/s]','interpreter','latex')
ylabel('$H_{sig}$ [m]','interpreter','latex')
title(['Waves on',' ',planet_to_run.name],'interpreter','latex');
grid on;


