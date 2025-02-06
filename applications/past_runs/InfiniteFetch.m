clc
clear
close all

% DEEPWATER WAVES FOR SEDIMENT ENTRAINMENT AT ONTARIO LACUS
addpath(fullfile('..','..','planetwaves'))  
addpath(fullfile('..','..','planetwaves','pre_analysis'))


% infinite fetch

% MODEL INPUTS
%lakes = {'Titan-CH3H8N2','Titan-OntarioLacus','Titan-LigeiaMare','Titan-CH4N2'};
lakes = {'Titan-CH3H8N2'};

zDep = 75.*ones(50,50);
buoy_loc = [25 25];                                                        % measure at deepest location
grid_resolution = [10*1000 10*1000];                                             % pixel width and pixel height [m]
test_speeds = 4;%[0.3:0.1:5];                                              % surface wind speeds [m/s]
time_to_run = 60*10;                                                          % time to run model [s]
wind_direction = 0;                                                        % wind direction [rad]

figure;
waveheight_ax = axes;
grid on;
xlabel('u10 [m/s]')
ylabel('sigH[m]')
hold on;
xlim([0 max(test_speeds)+0.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

    

    for i = 1:numel(test_speeds)


            Wind.speed = test_speeds(i);
            Model = calc_cutoff_freq(Planet,Model,Wind);

            [avgH, htgrid, ~, ~, ~,~, PeakWave] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
            plot(height_ax,1:numel(avgH),avgH,'-','DisplayName',num2str(Wind.speed))
            legend(height_ax,'show','Location','best')
            drawnow

            plot_height(c,i) = avgH(end);
            
  
            yline(height_ax,PeakWave.H(Model.long,Model.lat),'--k')
            
            H_mat{c,i} = PeakWave.H;
            H0(c,i) = PeakWave.H(Model.long,Model.lat);
            Cg0_mat{c,i} = PeakWave.cg;
            Cg0(c,i) = PeakWave.cg(Model.long,Model.lat);
            C0_mat{c,i} = PeakWave.c;
            C0(c,i) = PeakWave.c(Model.long,Model.lat);
            L_mat{c,i} = PeakWave.L;
            L0(c,i) = PeakWave.L(Model.long,Model.lat);
            T_mat{c,i} = PeakWave.T;
            T(c,i) = PeakWave.T(Model.long,Model.lat);

            %save(fullfile('..','past_run','Waves_03_43.mat')
            

    
         
    end
    


    
    pz = plot_height(c,:);
    plot(waveheight_ax,test_speeds(pz ~=0),pz(pz ~=0),'-s','LineWidth',2,'DisplayName',planet_to_run);
    
    legend(waveheight_ax,'show','Location','best')
    drawnow
    

        
end

figure; imagesc(PeakWave.H)