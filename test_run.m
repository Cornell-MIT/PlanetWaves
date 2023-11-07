clc
clear
close all

% % DATA SOURCE: https://www.ndbc.noaa.gov/station_realtime.php?station=45012
% LakeOntario45012 = readtable("LakeOntario_45012_5days.txt","TreatAsMissing","MM");
% LakeOntario45012 = renamevars(LakeOntario45012,['x_YY'],['YY']);
% 
% 
% pos_windspeed = unique(LakeOntario45012.WSPD);
% 
% 
% for ii = 1:length(pos_windspeed)
%     i = (LakeOntario45012.WSPD == pos_windspeed(ii));
%     a = LakeOntario45012(i,:);
%     avg_waveheight(ii) = mean(a.WVHT,"omitmissing");
%     std_waveheight(ii) = std(a.WVHT,"omitmissing");
% end
% 
% 
% LO_WVHT = dictionary(pos_windspeed,avg_waveheight');
% LO_WVHT_STD = dictionary(pos_windspeed,std_waveheight');

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% INPUT PARAMETERS:
% (1) PLANET CONDITIONS
%   (a) TITAN
Titan.rho_liquid = 540;                                                    % Hydrocarbon Liquid density [kg/m3]
Titan.nu_liquid = 3e-7;                                                    % Hydrocarbon Liquid Viscocity [m2/s]
Titan.nua = 0.0126/1e4;                                                    % Titan atmospheric gas viscocity [m2/s]
Titan.gravity = 1.352;                                                     % Titan Gravity [m/s2]
Titan.surface_temp = 92;                                                   % Titan Surface Temperature [K]
Titan.surface_press = 1.5*101300;                                          % Titan Surface Pressure [Pa]
Titan.surface_tension = 0.018;                                             % Hydrocarbon Liquid Surface Tension [N/m]
Titan.name = 'Titan';

%   (b) EARTH
Earth.rho_liquid = 997;                                                    % Water Liqid Density [kg/m3]         
Earth.nu_liquid = 1.143e-6;                                                % Water Liquid Viscocity [m2/s]
Earth.nua = 1.48e-5;                                                       % Earth atmospheric gas viscocity [m2/s]
Earth.gravity = 9.81;                                                      % Earth Gravity [m/s2]
Earth.surface_temp = 288;                                                  % Earth Surface Temperature [K]
Earth.surface_press = 1*101300;                                            % Earth Surface Pressure [Pa]
Earth.surface_tension = 0.072;                                             % Water Liquid Surface Tension [N/m]
Earth.name = 'Earth';

% (2a) MODEL GEOMETRY
Model.m = 20;                                                              % Number of Grid Cells in X-Dimension
Model.n = 100;                                                             % Number of Grid Cells in Y-Dimension
Model.o = 100;                                                             % Number of Frequency bins
Model.p = 288;                                                             % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.long = 10;                                                           % longitude grid point for sampling during plotting 
Model.lat = 10;                                                            % latitude grid point for sampling during plotting
Model.gridX = 10;                                                          % Grid size in X-dimension [m]
Model.gridY = 10;                                                          % Grid size in Y-dimension [m]
Model.mindelt = 0.0001;                                                    % minimum time step
Model.maxdelt = 2000.0;                                                    % maximum time step
Model.time_step = 10;                                                     % Maximum Size of time step [s]
Model.num_time_steps = 10;                                                 % Length of model run (in terms of # of time steps)
Model.tolH = NaN;                                                          % tolerance threshold for maturity 
Model.cutoff_freq = 50;                                                    % cutoff frequency bin from diagnostic to advection
Model.min_freq = 0.05;                                                     % minimum frequency to model
Model.max_freq = 35;                                                       % maximum frequency to model

% (2b) BUOY SPECIFC: 
% STATION 45012:
% https://www.ndbc.noaa.gov/station_page.php?station=45012
Model.z_data = 10;                                                        % elevation of wind measurement [m]
deep_bathy = 100.*ones(Model.m,Model.n);                                  % depth of water column beneath buoy [m]
Model.bathy_map = deep_bathy;                                             % Bathymetry of model basin [m]

% (2c) TUNING PARAMETERS 
Model.tune_A1 = 0.11;
Model.tune_mss_fac = 240;
Model.tune_Sdt_fac = 0.001;
Model.tune_Sbf_fac = 0.002;
Model.tune_cotharg = 0.2;
Model.tune_n = 2.4;

% (3) NEAR-SURFACE WIND CONDITIONS
Wind.speed = [3 5 7];                                                            % magnitude of incoming wind [m/s]
Wind.dir = 0;                                                              % direction of incoming wind [radians]

% (4) Unidirectional currents
Uniflow.East = 0;                                                          % eastward unidirectional current [m/s]
Uniflow.North = 0;                                                         % northward unidirectional current [m/s]

% (5) HOUSEKEEPING
Etc.showplots = 0;
Etc.savedata = 0;
Etc.showlog = 1;

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


planet_to_run = Earth;

[sigH,htgrid,freqspec] = makeWaves(planet_to_run,Model,Wind,Uniflow,Etc); 


UMWM_WVHT = dictionary(Wind.speed,sigH(:,end));

PM_H = 0.22.*(Wind.speed.^2)./(planet_to_run.gravity);

figure;
plot(Wind.speed,PM_H,'-b','LineWidth',3)
% plot(keys(LO_WVHT),values(LO_WVHT),'-ok')
% hold on
% plot(keys(LO_WVHT_STD),values(LO_WVHT) + values(LO_WVHT_STD),'--ok')
% plot(keys(LO_WVHT_STD),values(LO_WVHT) - values(LO_WVHT_STD),'--ok')
plot(keys(UMWM_WVHT),values(UMWM_WVHT),'-r')
%legend('Lake Ontario 45012 5-days','Lake Ontario + STD','Lake Ontario - STD','UMWM-Titan','location','best')
legend('Pierson-Moskowitz','UMWM-Titan')
grid on
xlabel('wind speed [m/s]')
ylabel('sig wave height [m]')
title(planet_to_run.name)

