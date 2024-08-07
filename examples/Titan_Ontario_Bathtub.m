clc
clear
close all

make_gif = 0;
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
% lakes = {'Titan-OntarioLacus','Titan-LigeiaMare','Titan-CH4N2','Titan-CH3H8N2'};
lakes = {'Titan-OntarioLacus'};
buoy_loc = [60 55];                                                        % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
test_speeds = [1:3];                                                     % wind speed
% time_to_run = 60*10;                                                        % time to run model
time_to_run = 60*10;
wind_direction = pi/2;                                                     % wind direction

figure;
waveheight_ax = axes;
grid on;
xlabel('u10 [m/s]')
ylabel('sig[m]')
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL MODEL
% degrade depth profile so model doesnt take as long to run
[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.2);
for composition = 1:numel(lakes)
    planet_to_run = lakes{composition};
    % populate model classes
    [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
    % update grid resolution
    Model.gridX = grid_resolution(1);                                              
    Model.gridY = grid_resolution(2);                                               
              
    Model.cutoff_freq = round((15/35)*Model.Fdim);
    
    figure;
    imagesc(zDep)
    hold on;
    plot(buoy_loc(1),buoy_loc(2),'or','MarkerFaceColor','r')
    colorbar;
    
    title('input bathymetry')
    [wx,wy] = pol2cart(Wind.dir,1);
    plot(Model.long,Model.lat,'pentagram','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20)
    quiver(Model.long,Model.lat, 5*wx, 5*wy, 'r', 'MaxHeadSize', 1);
    
    drawnow
    hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RUN MODEL
    
    
    lake_slope = 0.5e-3;
    d = Model.bathy_map(Model.long,Model.lat):-lake_slope:10^-4;
    % d50 = [6.35e-5 0.1];% diameters for [finegrain to 10 cm] [m]
    d50 = [6.35e-5];
    rho_s = 940; % ice grains
    
    
    figure;
    crash_ax = axes;
    xlabel('wind speed [m/s]')
    ylabel('entrainment depth [m]')
    hold on;
    title(planet_to_run)
    
    figure;
    height_ax = axes;
    xlabel('model time')
    ylabel('wave height')
    hold on;
    title(planet_to_run)
    
    for i = 1:numel(test_speeds)
        
    
            Wind.speed = test_speeds(i);
    
            [avgH, ~, Espec, ~, ~,~, PeakWave] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
            plot(height_ax,1:numel(avgH),avgH,'-','DisplayName',num2str(Wind.speed))
            legend('show','Location','best')
            drawnow

            plot_height(composition,i) = avgH(end);
            
    
            % shoal the waves 
            H0 = PeakWave.H(Model.long,Model.lat);
            Cg0 = PeakWave.cg(Model.long,Model.lat);
            C0 = PeakWave.c(Model.long,Model.lat);
            L0 = PeakWave.L(Model.long,Model.lat);
            T = PeakWave.T(Model.long,Model.lat);
            alpha_0 = 0; % assuming incoming wave crests parallel to shore contour

         for a = 1:numel(d50)
            for z = 1:numel(d)

                
                L_shoal(z) = L0*sqrt(tanh((4*(pi^2)*d(z)/(T^2*Planet.gravity))));
                d_L = d(z)/L_shoal(z);
                alpha(z) = asin(tanh((2*pi.*d_L).*sin(alpha_0)));
                C_shoal(z) = C0*tanh(2*pi*d_L);
                n(z) = 0.5*(1 + ((4*pi*d_L)/(sinh(4*pi*d_L))));
                Cg_shoal(z) = n(z)*C_shoal(z);
                KR(z) =  sqrt(cos(alpha_0)/cos(alpha(z)));
                KS(z) = sqrt(Cg0/Cg_shoal(z)); 
                H_shoal(z) = KR(z)*KS(z)*H0; 

                d0_shoal(z) = (H_shoal(z)/2).*((1)./sinh(2*pi*d_L)); % cosh @ z = -d = 1
                um_shoal(z) = (H_shoal(z)/2).*((Planet.gravity*T)/L_shoal(z)).*((1)./(cosh(2*pi*d_L))); % cosh @ z = -d = 1


                shields(z,i) = (Planet.rho_liquid.*(um_shoal(z).^2))./((rho_s - Planet.rho_liquid)*Planet.gravity.*d50(a));
                KM_crash(z,i) = 0.3*sqrt(d0_shoal(z)/d50(a));

                entrained = find(shields(:,i)>KM_crash(:,i),1,'first');

                if isempty(entrained) % does not reach entrainment
                    d_crash(i,a) = NaN;
                else % does reach entrainment
                    d_crash(i,a) = d(entrained);
                end

            end

    
            if a == 1
                figure;
                shoal_ax = axes;
                plot(shoal_ax,d,L_shoal,'--r','LineWidth',5)
                hold on
                plot(shoal_ax,d,C_shoal,'--g','LineWidth',5)
                plot(shoal_ax,d,n,'--b','LineWidth',5)
                plot(shoal_ax,d,Cg_shoal,'--m','LineWidth',5)
                plot(shoal_ax,d,H_shoal,'--c','LineWidth',5)
                legend('L','C','n','Cg','H','Location','best')
                title('shoaling')
                set(gca, 'XScale', 'log')
                set(gca,'Xdir','reverse')
                grid on;
                xlabel('depth [m]')
                title(['u = ' num2str(Wind.speed) ' m/s'])
                drawnow
            end



        end


        if d_crash(:,a) > 0
          plot(crash_ax,test_speeds,d_crash(:,a),'-s','LineWidth',1,'DisplayName',num2str(d50(a)))
          legend('show','Location','best')
          drawnow
        end
    end
        pz = plot_height(composition,:);
        plot(waveheight_ax,test_speeds(pz ~=0),pz(pz ~=0),'--s','LineWidth',2,'DisplayName',planet_to_run);
        legend('show','Location','best')
    
end
     


% [rows, columns, numberOfColorChannels] = size(zDep);
% 
% x_extent = 1:columns; 
% y_extent = 1:rows; 
% [X,Y] = meshgrid(x_extent,y_extent);
% 
% figure
% contour3(X,Y,zDep,[0:20:max(max(zDep))],'-k','ShowText','on','LabelSpacing',1000)
% hold on;
% [Mcon,Ccon] = contour3(X, Y, zDep,[d_crash(end,1) d_crash(end,1)],'Color',[62/255 150/255 81/255],'LineWidth',2,'ShowText',0); % entrainment depth for max wind speed for smallest grains
% [Mcon1,Ccon1] = contour3(X, Y, zDep,[d_crash(end,2) d_crash(end,2)],'Color',[204/255 37/255 41/255],'LineWidth',2,'ShowText',0); % entrainment depth for max wind speed for largest grains
% view(0,90)
% 
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% grid off;




