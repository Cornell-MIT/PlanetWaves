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
planet_to_run = 'Titan';
buoy_loc = [60 55];                                                        % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
test_speeds = 3;                                                     % wind speed
time_to_run = 60*2;                                                          % time to run model
wind_direction = pi/2;                                                        % wind direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL MODEL
% degrade depth profile so model doesnt take as long to run
[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.1);
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
hold off;
title('input bathymetry')
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL


lake_slope = 0.5e-3;
d = Model.bathy_map(Model.long,Model.lat):-lake_slope:10^-4;
d50 = [6.35e-5];% 0.1];% diameters for [finegrain to 10 cm] [m]
rho_s = 940; % ice grains


% Preallocate cell arrays to store results
% myHsig = cell(1, numel(test_speeds));
% htgrid = cell(1, numel(test_speeds));
% E_spec = cell(1, numel(test_speeds));
% Cg = cell(1,numel(test_speeds));
figure;
for a = 1:numel(d50)
    for i = 1:numel(test_speeds)
    
        Wind.speed = test_speeds(i);
        
        [~, htgrid{i}{a}, E_each{i}{a}, ~, ~,~, PeakWave] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
    
        
        % shoal the waves 
        H0 = PeakWave.H(Model.long,Model.lat);
        Cg0 = PeakWave.cg(Model.long,Model.lat);
        C0 = PeakWave.c(Model.long,Model.lat);
        L0 = PeakWave.L(Model.long,Model.lat);
        T = PeakWave.T(Model.long,Model.lat);
        alpha_0 = 0; % information in wavemodel? same as wind?
        
        
        for z = 1:numel(d)
        
            L_shoal(z) = L0*sqrt(tanh((4*(pi^2)*d(z)))/(Planet.gravity*(T^2))); % Eckart aprox 
            d_L = d(z)/L_shoal(z);
            alpha(z) = asin(tanh((2*pi.*d_L).*sin(alpha_0)));
            C_shoal(z) = C0*tanh(2*pi*d_L);
            n(z) = 0.5*(1 + ((4*pi*d_L)/(sinh(4*pi*d_L))));
            Cg_shoal(z) = n(z)*C_shoal(z);
            KR(z) =  sqrt(cos(alpha_0)/cos(alpha(z)));
            KS(z) = sqrt(Cg0/Cg_shoal(z)); 
            H_shoal(z) = KR(z)*KS(z)*H0; 
        
            d0_shoal(z) = H_shoal(z)./sinh(2*pi*d_L);
            um_shoal(z) = (pi*d0_shoal(z))/T;
        
            
            shields(z,i) = (Planet.rho_liquid.*(um_shoal(z).^2))./((rho_s - Planet.rho_liquid)*Planet.gravity.*d50(a));
            KM_crash(z,i) = 0.3*sqrt(d0_shoal(z)/d50(a));
        
            entrained = find(shields(:,i)>KM_crash(:,i),1,'first');
    
            if isempty(entrained) % does not reach entrainment
                d_crash(i) = NaN;
            else % does reach entrainment
                d_crash(i) = d(entrained);
            end
    
        end
    
    
    figure;
    plot(d,L_shoal,'--r','LineWidth',5)
    hold on
    plot(d,C_shoal,'--g','LineWidth',5)
    plot(d,n,'--b','LineWidth',5)
    plot(d,Cg_shoal,'--m','LineWidth',5)
    plot(d,H_shoal,'--c','LineWidth',5)
    legend('L','C','n','Cg','H')
    title('shoaling')
    set(gca, 'XScale', 'log')
    set(gca,'Xdir','reverse')
    grid on;
    xlabel('depth [m]')
    
    
    
    
    end
    
    

    plot(test_speeds,d_crash,'-s','LineWidth',1,'DisplayName',num2str(d50(a)))
    hold on
    xlabel('wind speed [m/s]')
    ylabel('entrainment depth [m]')
end
legend('show','Location','best')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  PLOT RESULTS
% 
% figure
% for i = 1:numel(test_speeds)
%     plot(myHsig{i}, '-', 'LineWidth', 3, 'DisplayName', num2str(test_speeds(i)))
%     hold on
% end
% grid on;
% legend('show', 'Location', 'northwest','interpreter','latex');
% title(['Waves on',' ',Planet.name],'interpreter','latex');
% xlabel('model time step [$\Delta$ t]','interpreter','latex')
% ylabel('significant wave height [m]','interpreter','latex')
% 
% 
% for k = 1:numel(test_speeds)
%     buoy_waves(k) = htgrid{k}{end}(Model.long,Model.lat);
% end
% 
% figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1])
% for speed = 1:numel(test_speeds)
% 
%     subplot(1,3,3)
%     plot(test_speeds,buoy_waves,'-sb','LineWidth',1,'MarkerFaceColor','b')
%     hold on
%     plot(test_speeds(speed),htgrid{speed}{end}(Model.long,Model.lat),'-sr','LineWidth',1,'MarkerFaceColor','r')
%     hold off;
%     xlabel('u [m/s]','interpreter','latex')
%     ylabel('Hsig [m]','interpreter','latex')
%     grid on;
% 
%     plot_grid = htgrid{speed}{end}';
%     plot_grid(isnan(plot_grid)) = 0;
%     plot_alpha_data = ones(size(plot_grid));
%     plot_alpha_data(plot_grid==0) = 0;
% 
%     ax1 = subplot(1,3,[1,2]);
%     h1 = imagesc(plot_grid);
%     colormap cool
%     xlabel('longitude [km]')
%     ylabel('latitude [km]')
%     title(sprintf('u = %i m/s',test_speeds(speed)))
%     c1 = colorbar;
%     c1.Label.String = 'Hsig [m]';
%     clim([0 1.2])
%     set(h1, 'AlphaData', plot_alpha_data);
%     hold on;
% 
%     %contour(Model.bathy_map,min(min(Model.bathy_map)):10:max(max(Model.bathy_map)),'-k','LineWidth',2)
%     contour(plot_grid,'-k','LineWidth',2)
% 
%     grid on
%     new_xtick = get(gca, 'XTick')*(Model.gridX)/1000;
%     new_ytick = get(gca, 'YTick')*(Model.gridY)/1000;
%     set(gca, 'XTick',  get(gca, 'XTick'), 'XTickLabel', arrayfun(@(x) sprintf('%d', x), new_xtick, 'UniformOutput', false));
%     set(gca, 'YTick',  get(gca, 'YTick'), 'YTickLabel', arrayfun(@(y) sprintf('%d', y), new_ytick, 'UniformOutput', false));
% 
% 
%     [wx,wy] = pol2cart(Wind.dir,1);
%     plot(Model.long,Model.lat,'pentagram','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20)
% 
%     for i = 1:Model.LonDim
%         for j = 1:Model.LatDim
%             if Model.bathy_map(j,i) <=0
%                 quiver(i, j, (speed/5)*wx, (speed/5)*wy, 'r', 'MaxHeadSize', 1);
%             end
%         end
%     end
%     set(ax1,'Ydir','reverse')
% 
%     drawnow
%     if speed == 1 
%         if make_gif
%             gif('Ontario_Titan_Bathtub.gif','DelayTime',1,'overwrite',true)
%         end
%     else
%         if make_gif
%             gif
%         end
%     end
% end
% 
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(test_speeds,buoy_waves,'-sb','LineWidth',1,'MarkerFaceColor','b')
% %ylim([0 3])
% xlabel('u [m/s]','Interpreter','latex')
% ylabel('Hsig [m]','Interpreter','tex')
% grid on;
% 
% % maturity = wave_age{1};
% % maturity(maturity>=0.83) = 1;
% % maturity(maturity<0.83) = 0;
% % 
% % figure
% % h2 = imagesc(maturity);
% % alpha_maturity = ones(size(Model.bathy_map));
% % alpha_maturity(Model.bathy_map <= 0) = 0;
% % set(h2, 'AlphaData', alpha_maturity);
% % 
% % cRange = caxis; 
% % hold on;
% % contour(Model.bathy_map,min(min(Model.bathy_map)):10:max(max(Model.bathy_map)),'-k','LineWidth',2)
% % caxis(cRange);
% % title('Maturity of Waves')
% 
% 
% 
% 
% 
% 
