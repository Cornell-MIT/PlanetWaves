% clc
% clear
% close all
% 
% 
lake_slope = 0.5e-3;
d50 = [6.35e-5 0.1];
rho_s = [800 940 1500]; % [organic-ice ice organic]
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

zDep_orig = zDep;
% MODEL INPUTS
lakes = {'Titan-CH3H8N2','Titan-OntarioLacus','Titan-LigeiaMare','Titan-CH4N2',};
buoy_loc = [60 55];                                                        % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
test_speeds = [0.5 1 2 3];                                                 % wind speed
time_to_run = 60*10;                                                     % time to run model
wind_direction = pi/2;                                                     % wind direction

figure;
waveheight_ax = axes;
grid on;
xlabel('u10 [m/s]')
ylabel('sigH[m]')
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL MODEL
% degrade depth profile so model doesnt take as long to run
[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.2);


for c = 1:numel(lakes)

    planet_to_run = lakes{c};
    % populate model classes
    [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
    % update grid resolution
    Model.gridX = grid_resolution(1);                                              
    Model.gridY = grid_resolution(2);                                               

    Model.cutoff_freq = round((15/35)*Model.Fdim);

    figure;
    height_ax = axes;
    xlabel('model time')
    ylabel('wave height')
    title(planet_to_run)
    hold on;

    if c == 1
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
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RUN MODEL FOR DEEPWATER WAVES

    d = Model.bathy_map(Model.long,Model.lat):-lake_slope:10^-4;

    for i = 1:numel(test_speeds)


            Wind.speed = test_speeds(i);

            [avgH, ~, Espec, ~, ~,~, PeakWave] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
            plot(height_ax,1:numel(avgH),avgH,'-','DisplayName',num2str(Wind.speed))
            legend(height_ax,'show','Location','best')
            drawnow

            plot_height(c,i) = avgH(end);


            % shoal the waves 
            H0(c,i) = PeakWave.H(Model.long,Model.lat);
            Cg0(c,i) = PeakWave.cg(Model.long,Model.lat);
            C0(c,i) = PeakWave.c(Model.long,Model.lat);
            L0(c,i) = PeakWave.L(Model.long,Model.lat);
            T(c,i) = PeakWave.T(Model.long,Model.lat);
            alpha_0 = 0; % assuming incoming wave crests parallel to shore contour

         for a = 1:numel(d50)

            for z = 1:numel(d)


                L_shoal(z) = L0(c,i)*sqrt(tanh((4*(pi^2)*d(z)/(T(c,i)^2*Planet.gravity))));
                d_L = d(z)/L_shoal(z);
                alpha(z) = asin(tanh((2*pi.*d_L).*sin(alpha_0)));
                C_shoal(z) = C0(c,i)*tanh(2*pi*d_L);
                n(z) = 0.5*(1 + ((4*pi*d_L)/(sinh(4*pi*d_L))));
                Cg_shoal(z) = n(z)*C_shoal(z);
                KR(z) =  sqrt(cos(alpha_0)/cos(alpha(z)));
                KS(z) = sqrt(Cg0(c,i)/Cg_shoal(z)); 
                H_shoal(z) = KR(z)*KS(z)*H0(c,i); 

                d0_shoal(z) = (H_shoal(z)/2).*((1)./sinh(2*pi*d_L)); % cosh @ z = -d = 1
                um_shoal(z) = (H_shoal(z)/2).*((Planet.gravity*T)/L_shoal(z)).*((1)./(cosh(2*pi*d_L))); % cosh @ z = -d = 1

                
                for s = 1:numel(rho_s)

                    shields(z,i) = (Planet.rho_liquid.*(um_shoal(z).^2))./((rho_s(s) - Planet.rho_liquid)*Planet.gravity.*d50(a));
                    KM_crash(z,i) = 0.3*sqrt(d0_shoal(z)/d50(a));
    
                    entrained = find(shields(:,i)>KM_crash(:,i),1,'first');
    
                    if isempty(entrained) % does not reach entrainment
                        d_crash{s,c}(i,a) = NaN;
                    else % does reach entrainment
                        d_crash{s,c}(i,a) = d(entrained);
                    end
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
                plot(shoal_ax,d,T.*ones(size(d)),'--k','LineWidth',5)
                legend('L','C','n','Cg','H','T','Location','best')
                title('shoaling')
                set(gca, 'XScale', 'log')
                set(gca,'Xdir','reverse')
                grid on;
                xlabel('depth [m]')
                title(['u = ' num2str(Wind.speed) ' m/s at ' planet_to_run])
                drawnow
            end



        end



    end

        pz = plot_height(c,:);
        plot(waveheight_ax,test_speeds(pz ~=0),pz(pz ~=0),'-s','LineWidth',2,'DisplayName',planet_to_run);
        legend(waveheight_ax,'show','Location','best')
        drawnow
end


icecolor = {'#42F2F7','#46ACC2','#498C8A','#4B6858'}; % ice colors (blues)
organicicecolor = {'#a53860','#da627d','#ffa5ab','#450920'}; % organic-ice (reds)
organiccolor =  {'#e9d700','#dab600','#a98600','#f8ed62'}; % organics (yellows)

figure
for c = 1:numel(lakes)
    for s = 1:numel(rho_s)
        if s == 1
            mycolor = organicicecolor;
        elseif s == 2
            mycolor = icecolor;
        elseif s == 3
            mycolor = organiccolor;
        end
        plot(test_speeds,d_crash{c}(:,1),'-o','Color',mycolor{c},'LineWidth',3,'MarkerFaceColor',mycolor{c},'DisplayName',[lakes{c} ' sand (rho = ' rho_s(s) ')'])
        hold on
        plot(test_speeds,d_crash{c}(:,2),'--s','Color',mycolor{c},'LineWidth',3,'MarkerFaceColor',mycolor{c},'DisplayName',[lakes{c} ' gravel(rho = ' rho_s(s) ')'])
    end
end
legend('show','Location','best')
grid on;
xlim([0 4])
xlabel('wind speed [m/s]')
ylabel('entrainment depth [m]')

% [rows, columns, numberOfColorChannels] = size(zDep_orig);
% 
% x_extent = 1:columns; 
% y_extent = 1:rows; 
% [X,Y] = meshgrid(x_extent,y_extent);
% 
% figure
% contour3(X,Y,zDep_orig,[1:10:max(max(zDep_orig))],'-k','ShowText','on','LabelSpacing',1000)
% hold on;
% [Mcon,Ccon] = contour3(X, Y, zDep_orig,[d_crash{2}(end,1) d_crash{2}(end,1)],'Color',[62/255 150/255 81/255],'LineWidth',2,'ShowText',0); % entrainment depth for max wind speed for smallest grains
% [Mcon1,Ccon1] = contour3(X, Y, zDep_orig,[d_crash{2}(end,2) d_crash{2}(end,2)],'Color',[204/255 37/255 41/255],'LineWidth',2,'ShowText',0); % entrainment depth for max wind speed for largest grains
% view(0,90)
% 
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% grid off;




% % c = 1;
% % figure
% % for i = [4 1 2 3]
% %     plot(test_speeds,d_crash{i}(:,1),'-o','Color',mycolor{c},'LineWidth',3,'MarkerFaceColor',mycolor{c},'DisplayName',[lakes{i} ' sand'])
% %     hold on
% %     plot(test_speeds,d_crash{i}(:,2),'--s','Color',mycolor{c},'LineWidth',3,'MarkerFaceColor',mycolor{c},'DisplayName',[lakes{i} ' gravel'])
% %     c = c + 1;
% % end
% % legend('show','Location','best')
% % grid on;
% % xlim([0 4])
% % xlabel('wind speed [m/s]')
% % ylabel('entrainment depth [m]')
% % c = 1;
% % for i = [4 1 2 3]
% %     pz = plot_height(i,:);
% %     plot(test_speeds(pz ~=0),pz(pz ~=0),'-s','Color',mycolor{c}, 'MarkerFaceColor',mycolor{c},'LineWidth',2,'DisplayName',lakes{i});
% %     hold on
% %     c = c + 1;
% % end
% % %legend('show','Location','best')
% [rows, columns, numberOfColorChannels] = size(zDep);
% 
% x_extent = 1:columns; 
% y_extent = 1:rows; 
% [X,Y] = meshgrid(x_extent,y_extent);
% 
% figure
% contour3(X,Y,zDep,[0.1:10:max(max(zDep))],'-k','ShowText','on','LabelSpacing',2000)
% hold on;
% [Mcon,Ccon] = contour3(X, Y, zDep,[d_crash{1}(end,1) d_crash{1}(end,1)],'Color',[62/255 150/255 81/255],'LineWidth',2,'ShowText',0); % entrainment depth for max wind speed for smallest grains
% [Mcon1,Ccon1] = contour3(X, Y, zDep,[d_crash{1}(end,2) d_crash{1}(end,2)],'Color',[204/255 37/255 41/255],'LineWidth',2,'ShowText',0); % entrainment depth for max wind speed for largest grains
% view(0,90)
% 
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% grid off;
