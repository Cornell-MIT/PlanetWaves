clc
clear
close all

% make plots of signifigant wave height for slopes

%% LOAD DATA:
% top row = methane alkane fraction (molar fraction of the alkanes (CH4+C2H6) that is methane)
% first column = temperature [K]
format long
density = load('../compositions/density.csv'); % [kg/m3]
kinematic_visc = load('../compositions/kinematic_visc.csv'); % [cm2/s]
dynamic_visc = load('../compositions/dynamic_visc.csv'); % [Pa*s]
surface_tension = load('../compositions/surface_tension.csv'); % [N/m]
methane_frac = load('../compositions/CH4_molefraction.csv'); % methane
ethane_frac = load('../compositions/C2H6_molefraction.csv'); % ethane
nitrogen_frac = load('../compositions/N2_molefraction.csv'); % nitrogen

%% INPUT PARAMETERS:
% (1) PLANET CONDITIONS
%   (a) TITAN
Titan.nua = 0.0126/1e4;                                                    % Titan atmospheric gas viscocity [m2/s]
Titan.gravity = 1.352;                                                     % Titan Gravity [m/s2]
Titan.surface_temp = 92;                                                   % Titan Surface Temperature [K]
Titan.surface_press = 1.5*101300;                                          % Titan Surface Pressure [Pa]
Titan.liquid.Hydrocarbon.rho_liquid = 540;                                 % Hydrocarbon liquid density [kg/m3]
Titan.liquid.Hydrocarbon.nu_liquid = 3e-7;                                 % Hydrocarbon liquid Viscocity [m2/s]
Titan.liquid.Hydrocarbon.sfct_liquid = 0.018;                              % Hydrocarbon liquid Surface Tension [N/m]

%   (a.1) VARYING METHANE:ETHANE:NITROGEN COMPOSITIONS @92 K               % Source: steckloff et al., 2020
temp_location = find(methane_frac(:,1)==Titan.surface_temp);               % Find row with temperature
methane = 0:0.25:1;                                                         % Fraction of methane
ethane = 1:-0.25:0;                                                         % Fraction of ethane
for i = 1:length(methane)
    nitrogen(i)=nitrogen_frac(temp_location,...                            % Fraction of nitrogen
        max(find(abs(nitrogen_frac(1,:) - methane(i))<0.001)));
    name(i) = string(sprintf('m%03.0fe%03.0fn%03.0f',...                 % Creating structural array with names
        100*methane(i),100*ethane(i),100*nitrogen(i)));
    Titan.liquid.(name(i)).rho_liquid = density(temp_location,...          % Liquid density [kg/m3]
        max(find(abs(density(1,:) - methane(i))<0.001))); 
    Titan.liquid.(name(i)).nu_liquid = kinematic_visc(temp_location,...    % Liquid viscocity [m2/s]
        max(find(abs(kinematic_visc(1,:) - methane(i))<0.001)))/10000;      
     Titan.liquid.(name(i)).sfct_liquid = surface_tension(temp_location,...% Liquid surface tension [N/m]
        max(find(abs(surface_tension(1,:) - methane(i))<0.001))); 
end
                                                      

%   (a.2) VARYING LAKE COMPOSITIONS @92 K
%   (a.2.1) TITAN: ONTARIO LACUS
%   51:38:11 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
Ontario_methane = 0.57;                                                    % Fraction of methane
Ontario_ethane = 0.43 ;                                                    % Fraction of ethane
Ontario_nitrogen = nitrogen_frac(temp_location,...                         % Fraction of nitrogen
        max(find(abs(nitrogen_frac(1,:) - Ontario_methane)<0.001)));
Titan.liquid.Ontario.rho_liquid = density(temp_location,...                % Ontario Lacus liquid density [kg/m3]
        max(find(abs(density(1,:) - Ontario_methane)<0.001))); 
Titan.liquid.Ontario.nu_liquid = kinematic_visc(temp_location,...          % Ontario Lacus liquid Viscocity [m2/s]
        max(find(abs(kinematic_visc(1,:) - Ontario_methane)<0.001)))/10000;      
Titan.liquid.Ontario.sfct_liquid = surface_tension(temp_location,...       % Ontario Lacus liquid Surface Tension [N/m]
        max(find(abs(surface_tension(1,:) - Ontario_methane)<0.001)));                                                       

%   (a.2.2) TITAN: PUNGA MARE
%   80:0:20 percent methane:ethane:nitrogen from Mastrogiuseppe 2018
Punga_methane = 1;                                                         % Fraction of methane
Punga_ethane = 0 ;                                                         % Fraction of ethane
Punga_nitrogen = nitrogen_frac(temp_location,...                           % Fraction of nitrogen
        max(find(abs(nitrogen_frac(1,:) - Punga_methane)<0.001)));
Titan.liquid.Punga.rho_liquid = density(temp_location,...                  % Punga Mare liquid density [kg/m3]
        max(find(abs(density(1,:) - Punga_methane)<0.001))); 
Titan.liquid.Punga.nu_liquid = kinematic_visc(temp_location,...            % Punga Mare liquid Viscocity [m2/s]
        max(find(abs(kinematic_visc(1,:) - Punga_methane)<0.001)))/10000;      
Titan.liquid.Punga.sfct_liquid = surface_tension(temp_location,...         % Punga Mare liquid Surface Tension [N/m]
        max(find(abs(surface_tension(1,:) - Punga_methane)<0.001)));                                                           

%   (a.2.3) TITAN: LIGIA MARE (MAX METHANE)
%   100:0:0 percent methane:ethane:nitrogen from LeGall 2016
LigiaMax_methane = 1;                                                      % Fraction of methane
LigiaMax_ethane = 0 ;                                                      % Fraction of ethane
LigiaMax_nitrogen = nitrogen_frac(temp_location,...                        % Fraction of nitrogen
        max(find(abs(nitrogen_frac(1,:) - LigiaMax_methane)<0.001)));
Titan.liquid.LigiaMax.rho_liquid = density(temp_location,...               % Ligia Mare (max) liquid density [kg/m3]
        max(find(abs(density(1,:) - LigiaMax_methane)<0.001))); 
Titan.liquid.LigiaMax.nu_liquid = kinematic_visc(temp_location,...         % Ligia Mare (max) liquid Viscocity [m2/s]
        max(find(abs(kinematic_visc(1,:) - LigiaMax_methane)<0.001)))/10000;      
Titan.liquid.LigiaMax.sfct_liquid = surface_tension(temp_location,...      % Ligia Mare (max) liquid Surface Tension [N/m]
        max(find(abs(surface_tension(1,:) - LigiaMax_methane)<0.001)));                                                            

%   (a.2.4) TITAN: LIGIA MARE (MIN METHANE)
%   50:39:10.73 percent methane:ethane:nitrogen from LeGall 2016
LigiaMin_methane = 0.56;                                                   % Fraction of methane
LigiaMin_ethane = 0.44 ;                                                   % Fraction of ethane
LigiaMin_nitrogen = nitrogen_frac(temp_location,...                        % Fraction of nitrogen
        max(find(abs(nitrogen_frac(1,:) - LigiaMin_methane)<0.001)));
Titan.liquid.LigiaMin.rho_liquid = density(temp_location,...               % Ligia Mare (min) liquid density [kg/m3]
        max(find(abs(density(1,:) - LigiaMin_methane)<0.001))); 
Titan.liquid.LigiaMin.nu_liquid = kinematic_visc(temp_location,...         % Ligia Mare (min) liquid Viscocity [m2/s]
        max(find(abs(kinematic_visc(1,:) - LigiaMin_methane)<0.001)))/10000;      
Titan.liquid.LigiaMin.sfct_liquid = surface_tension(temp_location,...      % Ligia Mare (min) liquid Surface Tension [N/m]
        max(find(abs(surface_tension(1,:) - LigiaMin_methane)<0.001)));                                                            
                                                  

%   (a.2.5) TITAN: KRAKEN MARE
%   70:16:14 percent methane:ethane:nitrogen from Poggiali et al., 2020
Kraken_methane = 0.81;                                                     % Fraction of methane
Kraken_ethane = 0.19 ;                                                     % Fraction of ethane
Kraken_nitrogen = nitrogen_frac(temp_location,...                          % Fraction of nitrogen
        max(find(abs(nitrogen_frac(1,:) - Kraken_methane)<0.001)));
Titan.liquid.Kraken.rho_liquid = density(temp_location,...                 % Kraken liquid density [kg/m3]
        max(find(abs(density(1,:) - Kraken_methane)<0.001))); 
Titan.liquid.Kraken.nu_liquid = kinematic_visc(temp_location,...           % Kraken liquid Viscocity [m2/s]
        max(find(abs(kinematic_visc(1,:) - Kraken_methane)<0.001)))/10000;      
Titan.liquid.Kraken.sfct_liquid = surface_tension(temp_location,...        % Kraken liquid Surface Tension [N/m]
        max(find(abs(surface_tension(1,:) - Kraken_methane)<0.001)));                                             

%   (b) EARTH
Earth.liquid.Water.rho_liquid = 997;                                       % Water liqid density [kg/m3]         
Earth.liquid.Water.nu_liquid = 1e-6;                                       % Water liquid viscocity [m2/s]
Earth.liquid.Water.sfct_liquid = 0.072;                                    % Water liquid surface ension [N/m]
Earth.gravity = 9.81;                                                      % Earth Gravity [m/s2]
Earth.surface_temp = 273;                                                  % Earth Surface Temperature [K]
Earth.surface_press = 1*101300;                                            % Earth Surface Pressure [Pa]
Earth.nua = 1.338/1e5;                                                     % Earth atmospheric gas viscocity [m2/s]

%% (2) MODEL GEOMETRY
Model.m = 30;                                                              % Number of Grid Cells in X-Dimension
Model.n = 30;                                                              % Number of Grid Cells in Y-Dimension
Model.o = 25;                                                              % Number of Frequency bins
Model.p = 288;                                                             % Number of angular (th) bins, must be factorable by 8 for octants to satisfy the Courant condition of numerical stability
Model.long = 15;                                                           % longitude grid point for sampling during plotting 
Model.lat = 15;                                                            % latitude grid point for sampling during plotting
Model.gridX = 1000.0;                                                      % Grid size in X-dimension [m]
Model.gridY = 1000.0;                                                      % Grid size in Y-dimension [m]
Model.mindelt = 0.0001;                                                    % minimum time step
Model.maxdelt = 2000.0;                                                    % maximum time step
Model.time_step = 100;                                                       % MAXIMUM SIZE OF TIME STEP [S]
Model.num_time_steps = 100;                                                 % LENGTH OF MODEL RUN (IN TERMS OF # OF TIME STEPS)
Model.tolH = NaN;                                                         % TOLERANCE THRESHOLD FOR MATURITY 

% % define a bathymetry with a constant slope
% bathy_map = ones(Model.m,Model.n);
% for i = 1:Model.m
%     bathy_map(i,:) = 50*(i/Model.m);
% end
% Model.bathy_map = bathy_map;                                             % Bathymetry of model basin [m]


% define a bathymetery map that is constant and deep
Model.bathy_map = 100.*ones(Model.m,Model.n);                              % Bathymetry of model basin [m]


%% (3) NEAR-SURFACE WIND CONDITIONS
Wind.dir = 0;                                                              % direction of incoming wind [radians]
%Wind.speed = log10(logspace(0,3,5));                                       % wind speed [m/s]
Wind.speed = [0 3];

%% (4) Unidirectional currents
Uniflow.East = 0;                                                          % eastward unidirectional current [m/s]
Uniflow.North = 0;                                                         % northward unidirectional current [m/s]

%% (5) HOUSEKEEPING
Etc.showplots = 0;
Etc.savedata = 0;
Etc.showlog = 0;

%% (6) RUN MODEL
Composition = string(fieldnames(Titan.liquid));

% Modeling Titan Waves
for liq = 1:length(Composition)
        [sigH{liq},htgrid{liq},E_each{liq}] = makeWaves_old(Titan,Titan.liquid.(Composition(liq)),Model,Wind,Uniflow,Etc);% [m]

        %Plotting individual liquid time step v sig H [TITAN]
        figure('WindowState','fullscreen')
        for i = 1:length(Wind.speed)
            plot(sigH{1,liq}(i,:),'DisplayName',['U = ',num2str(Wind.speed(1,i)),' m/s'],'LineWidth',.5)
            hold on
        end
        xlabel(['# time steps (delt = ', num2str(Model.time_step),' s)'])
        ylabel('sig H [m]')
        xlim([1 Model.num_time_steps])
        legend('Location','eastoutside','FontSize',12,'FontWeight','normal')
        title(['timestep vs sig H at ',num2str(Titan.surface_temp),' K for ', Composition{liq},' (methane:ethane:nitrogen)'],'FontSize',25)
        set(gca,'Fontsize',16,'FontWeight','bold')
        hold off
        saveas(gcf,[Composition{liq},'_timevSigH.png'])
        
        %Plotting individual liquid wind speed v sig H [TITAN]
        figure('WindowState','fullscreen') %plotting wind speed v sigH
        plot(Wind.speed,sigH{1,liq}(:,size(sigH{1,liq},2)),'DisplayName',Composition{liq},'LineWidth',1.5)
        legend('Location','eastoutside','FontSize',12,'FontWeight','normal')
        xlabel('windspeed [m/s]')
        ylabel('sig H [m]')
        xlim([0 max(Wind.speed)]);
        title(['windspeed vs sig H at ',num2str(Titan.surface_temp),' K for ', Composition{liq},' (methane:ethane:nitrogen)'],'FontSize',25)
        set(gca,'Fontsize',16,'FontWeight','bold')
        saveas(gcf,[Composition{liq},'_windvSigH.png'])
end

% Modeling Earth waves
[sigH_Water,htgrid_Water,E_each_Water] = makeWaves_old(Earth,Earth.liquid.Water,Model,Wind,Uniflow,Etc);% [m]

%Plotting individual water time step v sig H [EARTH]
figure('WindowState','fullscreen')
for i = 1:length(Wind.speed)
            plot(sigH_Water(i,:),'DisplayName',['U = ',num2str(Wind.speed(1,i)),' m/s'],'LineWidth',.5)
            hold on
end
xlabel(['# time steps (delt = ', num2str(Model.time_step),' s)'])
ylabel('sig H [m]')
xlim([1 Model.num_time_steps])
legend('Location','eastoutside','FontSize',12,'FontWeight','normal')
title(['timestep vs sig H at ',num2str(Earth.surface_temp),' K for Water'],'FontSize',25)
set(gca,'Fontsize',16,'FontWeight','bold')
hold off
saveas(gcf,'Water_timevSigH.png')

%Plotting individual water wind speed v sig H [EARTH]
figure('WindowState','fullscreen') %plotting wind speed v sigH
plot(Wind.speed,sigH_Water(:,size(sigH_Water,2)),'DisplayName','Water','LineWidth',1.5)
legend('Location','eastoutside','FontSize',12,'FontWeight','normal')
xlabel('windspeed [m/s]')
ylabel('sig H [m]')
xlim([0 max(Wind.speed)]);
title(['windspeed vs sig H at ',num2str(Earth.surface_temp),' K for Water'],'FontSize',25)
set(gca,'Fontsize',16,'FontWeight','bold')
saveas(gcf,'Water_windvSigH.png')
 
%Plotting varying methane:ethane:nitrogen wind speed v sig H + Earth
figure('WindowState','fullscreen')
for comp = 2:6
    plot(Wind.speed,sigH{1,comp}(:,size(sigH{1,comp},2)),'DisplayName', Composition{comp},'LineWidth',1.5)
    hold on
end
plot(Wind.speed,sigH_Water(:,size(sigH_Water,2)),'DisplayName', 'Water','LineWidth',1.5)
legend('Location','eastoutside','FontSize',12,'FontWeight','normal')
xlabel('windspeed [m/s]')
ylabel('sig H [m]')
xlim([0 max(Wind.speed)]);
title(['windspeed vs sig H at ',num2str(Titan.surface_temp), ' K for varying compositions (methane:ethane:nitrogen) (Water @ 288 K)'],'FontSize',25,'FontWeight','bold')
set(gca,'Fontsize',16,'FontWeight','bold')
saveas(gcf,'VaryingComps_windvSigH.png')

%Plotting varying lakes wind speed v sig H + Earth
figure('WindowState','fullscreen')
for lake = 7:11
    plot(Wind.speed,sigH{1,lake}(:,size(sigH{1,lake},2)),'DisplayName', Composition{lake},'LineWidth',1.5)
    hold on
end
%plot(Wind.speed,sigH_Water(:,size(sigH_Water,2)),'DisplayName', 'Water','LineWidth',1.5)
legend('Location','eastoutside','FontSize',12,'FontWeight','normal')
xlabel('windspeed [m/s]')
ylabel('sig H [m]')
xlim([0 max(Wind.speed)]);
title(['windspeed vs sig H at ',num2str(Titan.surface_temp), ' K for varying lakes (Water @ 288K)'],'FontSize',25,'FontWeight','bold')
set(gca,'Fontsize',16,'FontWeight','bold')
saveas(gcf,'TitanLakes_windvSigH.png')

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





