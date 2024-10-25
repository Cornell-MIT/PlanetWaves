clc
clear
close all

% DEEPWATER WAVES FOR SEDIMENT ENTRAINMENT
addpath(fullfile('..','planetwaves'))  
addpath(fullfile('..','planetwaves','pre_analysis'))
addpath(fullfile('..','data','Titan','TitanLakes','Bathymetries','bathtub_bathy'))
load('..\data\Titan\TitanLakes\Bathymetries\bathtub_bathy\ol_bathtub_0.002000_slope','zDep');

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


buoy_loc = [577 835];                                                      % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
test_speeds = [0.3:0.1:4.5];
time_to_run = 10*60;                                                          % time to run model
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
            Model = calc_cutoff_freq(Planet,Model,Wind);

            [avgH, ~, ~, ~, ~,~, PeakWave] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
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

            %save('Waves_03_43.mat')
    
         
    
         end
    


    
    pz = plot_height(c,:);
    plot(waveheight_ax,test_speeds(pz ~=0),pz(pz ~=0),'-s','LineWidth',2,'DisplayName',planet_to_run);
    legend(waveheight_ax,'show','Location','best')
    drawnow
    

        
end

