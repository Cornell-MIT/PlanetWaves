clc
clear
close all

addpath(fullfile('..','planetwaves'))  
addpath(fullfile('..','planetwaves','pre_analysis'))
addpath(fullfile('..','data','Titan','TitanLakes','Bathymetries','bathtub_bathy'))
load('..\data\Titan\TitanLakes\Bathymetries\bathtub_bathy\ol_bathtub_0.002000_slope','zDep');

% sediment entrainment in ontario lacus using planetwave

lakecolors = {'#F2D1C9','#E086D3','#8332AC','#462749'};
max_steepness = 1/7;
min_depth = 10^-4;
lake_slope = 0.002000;
d50 = [6.35e-5:1e-5:0.1];
rho_s = [800 940 1500]; % [organic-ice ice organic]
zDep_orig = zDep;
% MODEL INPUTS
lakes = {'Titan-CH3H8N2','Titan-OntarioLacus','Titan-LigeiaMare','Titan-CH4N2'};
ol = 2; lm = 3;
org = 1; ice = 2; fluffy = 3;

load('TitanLakesWaves.mat','test_speeds')
buoy_loc = [577 835];                                                      % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
%test_speeds = [0.3:0.1:4.5];
time_to_run = 60*10;                                                          % time to run model
wind_direction = 0;                                                        % wind direction

figure;
waveheight_ax = axes;
grid on;
xlabel('u10 [m/s]')
ylabel('sigH[m]')
hold on;
xlim([0 3.5])

figure;
shoal_ax = axes;
grid on;
xlabel('depth')
ylabel('wave height')
hold on;
set(shoal_ax, 'XDir','reverse')
set(shoal_ax,'XScale','log')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL MODEL
% degrade depth profile so model doesnt take as long to run
[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.03);

zDep(:,end) = [];

d_crash = {};
for c = 1:numel(lakes)

    planet_to_run = lakes{c};
    % populate model classes
    [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
    % update grid resolution
    Model.gridX = grid_resolution(1);                                              
    Model.gridY = grid_resolution(2);    
    

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

    d = max(max(Model.bathy_map)):-lake_slope:min_depth;

    for i = 1:numel(test_speeds)


            Wind.speed = test_speeds(i);
            % Model = calc_cutoff_freq(Planet,Model,Wind);

            % [avgH, ~, ~, ~, ~,~, PeakWave] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
            % plot(height_ax,1:numel(avgH),avgH,'-','DisplayName',num2str(Wind.speed))
            % legend(height_ax,'show','Location','best')
            % drawnow
            % 
            % plot_height(c,i) = avgH(end);
            % 
            % % shoal the waves 
            % H0(c,i) = PeakWave.H(Model.long,Model.lat);
            % Cg0(c,i) = PeakWave.cg(Model.long,Model.lat);
            % C0(c,i) = PeakWave.c(Model.long,Model.lat);
            % L0(c,i) = PeakWave.L(Model.long,Model.lat);
            % T(c,i) = PeakWave.T(Model.long,Model.lat);
            % save('LongRun.mat','C0','Cg0','H0','L0','plot_height','T','test_speeds','time_to_run')
            load('TitanLakesWaves.mat')

            alpha_0 = 0; % assuming incoming wave crests parallel to shore contour



         if H0(c,i) > 0

             
    
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
                    um_shoal(z) = (H_shoal(z)/2).*((Planet.gravity*T(c,i))/L_shoal(z)).*((1)./(cosh(2*pi*d_L))); % cosh @ z = -d = 1
                    break_frac(z) = H_shoal(z)/L_shoal(z);
                    break_frac(break_frac>=max_steepness) = NaN;
                end
                   
                
                    plot(shoal_ax,d(~isnan(break_frac)),H_shoal(~isnan(break_frac)),'LineWidth',3,'Color',lakecolors{c},'DisplayName',num2str(Wind.speed))
                    hold on;
                
                for s = 1:numel(rho_s)
    
                        for a = 1:numel(d50)
             
                            shields(:,i) = (Planet.rho_liquid.*(um_shoal.^2))./((rho_s(s) - Planet.rho_liquid)*Planet.gravity.*d50(a));
                            KM_crash(:,i) = 0.3*sqrt(d0_shoal/d50(a));
            
                            % figure
                            % plot(d,shields(:,i),'-k')
                            % hold on
                            % plot(d,KM_crash(:,i),'-r')
                            entrained_depth_index = find(shields(:,i)>KM_crash(:,i),1,'first');
    
                            % if ~isempty(entrained_depth_index)
                            %     xline(d(entrained_depth_index))
                            % end
    
                            if isempty(entrained_depth_index) % does not reach entrainment
                                d_crash{s,c}(i,a) = NaN;
                            else % does reach entrainment
                                if H_shoal(entrained_depth_index)/L_shoal(entrained_depth_index) < max_steepness
                                    d_crash{s,c}(i,a) = d(entrained_depth_index);
                                else % does not reach entrainment before wave breaking
                                    d_crash{s,c}(i,a) = -1;
                                end
                            end
                        end
                 end
    
                % if a == 1 &&  exist('L_shoal','var')
                %     figure;
                %     shoal_ax = axes;
                %     plot(shoal_ax,d,L_shoal,'-r','LineWidth',5)
                %     hold on
                %     plot(shoal_ax,d,C_shoal,'-g','LineWidth',5)
                %     plot(shoal_ax,d,n,'-b','LineWidth',5)
                %     plot(shoal_ax,d,Cg_shoal,'-m','LineWidth',5)
                %     plot(shoal_ax,d,H_shoal,'-c','LineWidth',5)
                %     plot(shoal_ax,d,T(c,i).*ones(size(d)),'-k','LineWidth',5)
                %     legend('L','C','n','Cg','H','T','Location','best')
                %     title('shoaling')
                %     set(gca, 'XScale', 'log')
                %     set(gca,'Xdir','reverse')
                %     grid on;
                %     xlabel('depth [m]')
                %     title(['u = ' num2str(Wind.speed) ' m/s at ' planet_to_run])
                %     drawnow
                % end
    
         
    
         end
    


    end
    % pz = plot_height(c,:);
    % plot(waveheight_ax,test_speeds(pz ~=0),pz(pz ~=0),'-s','LineWidth',2,'DisplayName',planet_to_run);
    % legend(waveheight_ax,'show','Location','best')
    % drawnow
    

        
end

figure;
for c = 1:numel(lakes)
    plot(d50,d_crash{ice,c}(end,:),'-','LineWidth',3,'Color',lakecolors{c},'DisplayName',lakes{c})
    hold on
end
legend('show')
xlabel('d50')
ylabel('entrainment depth')
hold off

icecolor = {'#42F2F7','#46ACC2','#498C8A','#4B6858'}; % ice colors (blues)
organicicecolor = {'#ffa5ab','#da627d','#a53860','#450920'}; % organic-ice (reds)
organiccolor =  {'#f8ed62','#e9d700','#dab600','#a98600'}; % organics (yellows)

figure
for s = 1:numel(rho_s)
    
    if s == 1
        mycolor = organicicecolor;
    elseif s == 2
        mycolor = icecolor;
    elseif s == 3
        mycolor = organiccolor;
    end

    for c = 1:numel(lakes)
  
        depth_entrain_sand = d_crash{s,c}(:,1);
        depth_entrain_sand(depth_entrain_sand<=0) = NaN;
        depth_entrain_grav = d_crash{s,c}(:,end);
        depth_entrain_grav(depth_entrain_grav<=0) = NaN;
        plot(test_speeds,depth_entrain_sand,'-^','Color',mycolor{c},'LineWidth',3,'MarkerFaceColor',mycolor{c},'DisplayName',[lakes{c} ' sand (rho = ' num2str(rho_s(s)) ')'])
        hold on
        plot(test_speeds,depth_entrain_grav,'-s','Color',mycolor{c},'LineWidth',3,'MarkerFaceColor',mycolor{c},'DisplayName',[lakes{c} ' gravel(rho = ' num2str(rho_s(s)) ')'])
    end

    % legend('show','Location','best')
    grid on;
    xlim([0 3.5])
    xlabel('wind speed [m/s]')
    ylabel('entrainment depth [m]')

end


[rows, columns, numberOfColorChannels] = size(zDep_orig);

x_extent = 1:columns; 
y_extent = 1:rows; 
[X,Y] = meshgrid(x_extent,y_extent);

figure
contour3(X,Y,zDep_orig,0:max(max(zDep_orig))/10:max(max(zDep_orig)),'-k','ShowText','on','LabelSpacing',1000);
hold on;
[sand_ol,~] = contour3(X, Y, zDep_orig,[d_crash{ice,ol}(end,1) d_crash{ice,ol}(end,1)],'Color',[62/255 150/255 81/255],'LineWidth',2,'ShowText',0); % entrainment depth for max wind speed for smallest grains
[grav_ol,~] = contour3(X, Y, zDep_orig,[d_crash{ice,ol}(end,end) d_crash{ice,ol}(end,end)],'Color',[204/255 37/255 41/255],'LineWidth',2,'ShowText',0); % entrainment depth for max wind speed for largest grains
view(0,90)

set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
grid off;



load('..\data\Titan\TitanLakes\Bathymetries\bathtub_bathy\lm_bathtub_0.002000_slope','zDep');

zDep_orig2 = zDep;
[rows, columns, numberOfColorChannels] = size(zDep_orig2);

x_extent = 1:columns; 
y_extent = 1:rows; 
[X,Y] = meshgrid(x_extent,y_extent);

figure
contour3(X,Y,zDep_orig2,0:max(max(zDep_orig2))/5:max(max(zDep_orig2)),'-k','ShowText','on','LabelSpacing',1000)
hold on;
[sand_lm,~] = contour3(X, Y, zDep_orig2,[d_crash{ice,lm}(end,1) d_crash{ice,lm}(end,1)],'Color',[62/255 150/255 81/255],'LineWidth',2,'ShowText',0); % entrainment depth for max wind speed for smallest grains
[grav_lm,~] = contour3(X, Y, zDep_orig2,[d_crash{ice,lm}(end,end) d_crash{ice,lm}(end,end)],'Color',[204/255 37/255 41/255],'LineWidth',2,'ShowText',0); % entrainment depth for max wind speed for largest grains
view(0,90)

set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
grid off;

sand_ol_contour = getContourLineCoordinates(sand_ol);
grav_ol_contour = getContourLineCoordinates(grav_ol);
sand_lm_contour = getContourLineCoordinates(sand_lm);
grav_lm_contour = getContourLineCoordinates(grav_lm);



load('mag_ang.mat')

[wind_speeds, counts] = unique_counts(round(mag_wind * 2) / 2);

figure
plot(wind_speeds,counts,'LineWidth',3)

counts(wind_speeds<0.5) = 0;

hold on
plot(wind_speeds,counts,'LineWidth',3)

xlabel('wind speed')
ylabel('count')

figure;
for  c = 1:numel(lakes)

    for  i = 1:numel(test_speeds)
        deepest_entrain(i) = d_crash{ice,c}(i,1)/d_crash{ice,c}(end,1);
        if i == 1
            wolman_miller(i) = deepest_entrain(i).*counts(2);
        elseif i == 2
            wolman_miller(i) = deepest_entrain(i).*counts(3);
        elseif i == 3
            wolman_miller(i) = deepest_entrain(i).*counts(5);
        elseif i == 4
            wolman_miller(i) = deepest_entrain(i).*counts(7);
        end
    end


    plot([0 test_speeds],[0 wolman_miller],'Color',lakecolors{c},'LineWidth',3,'DisplayName',lakes{c})
    hold on
end
xlim([0 3.5])
ylabel('freq x depth')
xlabel('wind speed')

function [uniqueVals, counts] = unique_counts(A)
    % Find the unique values in the array A
    uniqueVals = unique(A);
    
    % Count occurrences of each unique value
    counts = histc(A, uniqueVals); 
end