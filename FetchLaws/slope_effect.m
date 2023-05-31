clc
clear
close all

% make plots of signifigant wave height for slopes

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

%   (b) EARTH
Earth.rho_liquid = 997;                                                    % Water Liqid Density [kg/m3]         
Earth.nu_liquid = 1e-6;                                                    % Water Liquid Viscocity [m2/s]
Earth.gravity = 9.81;                                                      % Earth Gravity [m/s2]
Earth.surface_temp = 273;                                                  % Earth Surface Temperature [K]
Earth.surface_press = 1*101300;                                            % Earth Surface Pressure [Pa]
Earth.surface_tension = 0.072;                                             % Water Liquid Surface Tension [N/m]

% (2) MODEL GEOMETRY
Model.m = 30;                                                              % Number of Grid Cells in X-Dimension
Model.n = 30;                                                              % Number of Grid Cells in Y-Dimension
Model.o = 25;                                                              % Number of Frequency bins
Model.p = 288;                                                             % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.long = 15;                                                           % longitude grid point for sampling during plotting 
Model.lat = 15;                                                             % latitude grid point for sampling during plotting
Model.gridX = 1000.0;                                                      % Grid size in X-dimension [m]
Model.gridY = 1000.0;                                                      % Grid size in Y-dimension [m]
Model.mindelt = 0.0001;                                                    % minimum time step
Model.maxdelt = 2000.0;                                                    % maximum time step
Model.time_step = 1000;                                                       % Maximum Size of time step [s]
Model.num_time_steps = 10;                                               % Length of model run (in terms of # of time steps)
Model.tolH = NaN;                                                          % tolerance threshold for maturity 

% define a bathymetry with a constant slope
bathy_map = ones(Model.m,Model.n);
for i = 1:Model.m
    bathy_map(i,:) = 50*(i/Model.m);
end

Model.bathy_map = bathy_map;                                               % Bathymetry of model basin [m]

% (3) NEAR-SURFACE WIND CONDITIONS
Wind.dir = 0;                                                              % direction of incoming wind [radians]
Wind.speed = [1 5];                                                        % magnitude of incoming wind [m/s]

% (4) Unidirectional currents
Uniflow.East = 0;                                                          % eastward unidirectional current [m/s]
Uniflow.North = 0;                                                         % northward unidirectional current [m/s]

% (5) HOUSEKEEPING
Etc.showplots = 0;
Etc.savedata = 0;
Etc.showlog = 1;

% RUN MODEL
[sigH,htgrid,E_each] = makeWaves(Titan,Model,Wind,Uniflow,Etc); % [m]

figure;
plot(sigH(1,:))
hold on
plot(sigH(2,:))
plot(sigH(3,:))
plot(sigH(4,:))
plot(sigH(5,:))
plot(sigH(6,:))
plot(sigH(7,:))
plot(sigH(8,:))
plot(sigH(9,:))
plot(sigH(10,:))
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
% older runs:

% p = 64;
% o = 25;
% dr = pi/180;                                                               % conversion from degrees to radians
% dthd = 360/(p);
% th = rad2deg(([0:p-1]-p/2+0.5)*dthd*dr);
% f1 = 0.05;                                                                 % minimum frequency
% f2 = 35;                                                                   % maximum frequency
% % create frequency limits for spectrum
% dlnf=(log(f2)-log(f1))/(o-1);                                              % frequency step size for log normal distribution
% f = exp(log(f1)+[0:o-1]*dlnf);                                             % frequencies for spectrum
% dom = 2*pi*dlnf.*f;   
% 
% for q = 1:2
%     for i = 1:30
%         for j = 1:30
%             KK = squeeze(E{1,100}(i,j,:,:));
%             [max_num, max_idx]=max(KK(:));
%             [xxx,~] = ind2sub(size(KK),find(KK==max_num,1,'first'));
%             dom_w{q}(i,j) = dom(xxx);
%         end
%     end
% end

% 
% clc
% clear
% close all
% 
% % make plots of signifigant wave height for slopes
% 
% windspeeds = [0.54 1 3.3];
% 
% m = 50;                                                                    % number of gridpoints in x-direction
% n = 50;                                                                    % number of gridpoints in y-direction
% 
% 
% bathy_map = ones(m,n);
% for i = 1:m
%     bathy_map(i,:) = 50*(i/m);
% end
% 
% % surf(bathy_map)
% 
% 
% 
% rho_liquid = 590;
% nu_liquid = 7.5e-7; % m2/s
% planet_gravity = 1.352;
% planet_temp = 92;
% planet_press = 1.5*101300;
% surface_tension = 0.018;
% gridX = 1000.0;
% gridY = 1000.0;
% time_step = 1;
% num_time_steps = 1000;
% wind_dir = 0;
% 
% 
% [sigH,htgrid,E_each] = makeWaves(windspeeds,wind_dir,rho_liquid,nu_liquid,planet_gravity,planet_temp,planet_press,surface_tension,bathy_map,gridX,gridY,time_step,num_time_steps); % [m]
% 
% % hh_profile{1} = htgrid{1}(1:end,5:25);
% % hh_profile{2} = htgrid{2}(1:end,5:25);
% % hh_profile{3} = htgrid{3}(1:end,5:25);
% % hh_profile{4} = htgrid{4}(1:end,5:25);
% % 
% % subplot(2,1,1)
% % plot(1:length(bathy_map),-bathy_map)
% % hold off;
% % title('bathy map')
% % xlabel('distance from shore')
% % subplot(2,1,2)
% % plot(hh_profile{1}(:,10))
% % hold on;
% % plot(hh_profile{2}(:,10))
% % plot(hh_profile{3}(:,10))
% % plot(hh_profile{4}(:,10))
% % legend('1','3','5','10')
% % view(90,0)
% 
% % p = 64;
% % o = 25;
% % dr = pi/180;                                                               % conversion from degrees to radians
% % dthd = 360/(p);
% % th = rad2deg(([0:p-1]-p/2+0.5)*dthd*dr);
% % f1 = 0.05;                                                                 % minimum frequency
% % f2 = 35;                                                                   % maximum frequency
% % % create frequency limits for spectrum
% % dlnf=(log(f2)-log(f1))/(o-1);                                              % frequency step size for log normal distribution
% % f = exp(log(f1)+[0:o-1]*dlnf);                                             % frequencies for spectrum
% % dom = 2*pi*dlnf.*f;   
% % 
% % for q = 1:2
% %     for i = 1:30
% %         for j = 1:30
% %             KK = squeeze(E{1,100}(i,j,:,:));
% %             [max_num, max_idx]=max(KK(:));
% %             [xxx,~] = ind2sub(size(KK),find(KK==max_num,1,'first'));
% %             dom_w{q}(i,j) = dom(xxx);
% %         end
% %     end
% % end
% 
% % 
% 
% 
% 
% clc
% clear
% close all
% 
% % Lake Ontario Bathymetry: https://www.ncei.noaa.gov/products/great-lakes-bathymetry
% % Citation: National Geophysical Data Center, 1999. Bathymetry of Lake Ontario. National Geophyiscal Data Center, NOAA. https://doi.org/10.7289/V56H4FBH
% addpath('..\bathymetries\')
% % export DEM values into array for UMWM bathymetry
% D = imread('OntarioDEM1.tif');
% D(D==D(1,1)) = -1; % set land values to 0 (depth is positive)
% % figure;
% % contourf(D,20)
% % xlabel('x')
% % ylabel('y')
% % colormap(hot)
% % xlim([0 249])
% % ylim([0 70])
% 
% 
% dd = D(45:55,190:200);
% contourf(dd)
% 
% [m,n] = size(dd);
% 
% windspeeds = 7.7; % 15 knots (waves build to 0.6-1.2 meters)
% 
% rho_liquid = 997;
% nu_liquid = 1e-6;
% planet_gravity = 9.81;
% planet_temp = 273;
% planet_press = 1*101300;
% surface_tension =  0.072;
% bathy_map = dd;
% gridX = 1000.0;
% gridY = 1000.0;
% time_step = 1;
% num_time_steps = 100;
% 
% 
% wind_dir = 0;
% 
% 
% [sigH,htgrid,E] = makeWaves(windspeeds,wind_dir,rho_liquid,nu_liquid,planet_gravity,planet_temp,planet_press,surface_tension,bathy_map,gridX,gridY,time_step,num_time_steps); % [m]

% clc
% clear
% close all
% 
% % contour plots of wavenumber vs angle 
% addpath 'C:\Users\schne\OneDrive\Desktop\Main\Work\MIT\Github_Repos\umwm_titan\FetchLaws\Titan'
% p = 64;
% o = 25;
% dr = pi/180;                                                               % conversion from degrees to radians
% dthd = 360/(p);
% th = rad2deg(([0:p-1]-p/2+0.5)*dthd*dr);
% f1 = 0.05;                                                                 % minimum frequency
% f2 = 35;                                                                   % maximum frequency
% % create frequency limits for spectrum
% dlnf=(log(f2)-log(f1))/(o-1);                                              % frequency step size for log normal distribution
% f = exp(log(f1)+[0:o-1]*dlnf);                                             % frequencies for spectrum
% dom = 2*pi*dlnf.*f;                                                        % angular frequency (w = 2pi*f)
% idx = 1;
% 
% [rr,tt] = meshgrid(th,dom);
% for ii = 1:100
%     load(['Titan_1_' num2str(ii) '.mat'])
%     contourf(th,dom,squeeze(E(5,5,:,:)))
%     colorbar
%     set(gca, 'YScale', 'log')
%     xlabel('direction [deg]')
%     ylabel('angular frequency [rad/s]')
%     title(['E Time Step: ' num2str(ii)])
%     ylim([dom(1) dom(end)])
%     xlim([th(1) th(end)])
% %     xx = rr.*cos(tt);
% %     yy = rr.*sin(tt);
% %     hh = polar(xx,yy);
% %     hold on
% %     contour(xx,yy,squeeze(E(5,5,:,:)))
% %     set(hh,'Visible','off')
% %     axis off
% %     axis image
%     drawnow
% %     frame = getframe(gcf);
% %     im{idx} = frame2im(frame);
% %     idx = idx + 1;
% end
% % filename = 'E.gif';
% % idx_end = idx;
% % for idx = 1:idx_end-1
% %     [A,map] = rgb2ind(im{idx},256);
% % if idx == 1
% %     imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
% % else
% %     imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
% % end
% %     pause(0.1)
% % end





