clc
clear
close all

% DEEPWATER WAVES FOR SEDIMENT ENTRAINMENT AT ONTARIO LACUS
addpath(genpath(fullfile('..','..','planetwaves')))
addpath(fullfile('..','..','data','Titan','TitanLakes','Bathymetries','bathtub_bathy'))
load('..\..\data\Titan\TitanLakes\Bathymetries\bathtub_bathy\ol_bathtub_0.002000_slope','zDep');

save_file = 'past_runs/Titan_DeepwaterWaves.mat'; % name of file with saved data
if ~exist('past_runs', 'dir')
    mkdir('past_runs')
else
    disp('deleted old run')
    delete('./past_runs/Titan_DeepwaterWaves.mat')
end

% MODEL INPUTS
lakes = {'Titan-CH3H8N2','Titan-OntarioLacus','Titan-LigeiaMare','Titan-CH4N2'};

[~,I] = max(zDep,[],"all","linear");
[by, bx] = ind2sub(size(zDep),I);
buoy_loc = [bx by];                                                        % measure at deepest location
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
test_speeds = [0.3:0.1:4.6];                                               % surface wind speeds [m/s]
time_to_run = 60*10;                                                       % time to run model [s]
wind_direction = deg2rad(-45);                                             % wind direction [rad], towards annuli in Barnes+2014

figure("Name",'PlanetWaves: u vs H');
waveheight_ax = axes;
grid on;
xlabel('u10 [m/s]')
ylabel('sigH[m]')
hold on;
xlim([0 max(test_speeds)+0.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL MODEL
% degrade depth profile so model doesnt take as long to run
[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.03);
% remove empty column
zDep(:,end) = [];



for c = 1:numel(lakes)

    planet_to_run = lakes{c};
    % populate model classes
    [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);

    % update grid resolution
    Model.gridX = grid_resolution(1);                                              
    Model.gridY = grid_resolution(2);    

    Model.Dirdim = 8;


    if c == 1
        make_input_map(Planet,Model,Wind)
    end

    figure('Name',['Growth in ',planet_to_run]);
    height_ax = axes;
    xlabel('model time (min)')
    ylabel('wave height (m)')
    title(planet_to_run)
    hold on;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RUN MODEL FOR DEEPWATER WAVES

    for i = 1:numel(test_speeds)


            Wind.speed = test_speeds(i);
            Model = calc_cutoff_freq(Planet,Model,Wind);

            [avgH, htgrid, ~, ~, ~,~, PeakWave] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  

            if sum(avgH) ~= 0
                plot(height_ax,1:numel(avgH),avgH,'-','DisplayName',num2str(Wind.speed))
                legend(height_ax,'show','Location','best')
                drawnow
                plot_height(c,i) = avgH(end);
                % take value at tallest wave location
                [H_mat{c,i},ii] =  max(PeakWave.H,[],"all","omitnan","linear");
                H0(c,i) = PeakWave.H(ii);
                Cg0_mat{c,i} = PeakWave.cg;
                Cg0(c,i) = PeakWave.cg(ii);
                C0_mat{c,i} = PeakWave.c;
                C0(c,i) = PeakWave.c(ii);
                L_mat{c,i} = PeakWave.L;
                L0(c,i) = PeakWave.L(ii);
                T_mat{c,i} = PeakWave.T;
                T(c,i) = PeakWave.T(ii);
            else
                H_mat{c,i} = 0;
                H0(c,i) = 0;
                Cg0_mat{c,i} = 0;
                Cg0(c,i) = 0;
                C0_mat{c,i} = 0;
                C0(c,i) = 0;
                L_mat{c,i} = 0;
                L0(c,i) =0;
                T_mat{c,i} = NaN;
                T(c,i) = NaN;
            end
            save(save_file,'test_speeds','time_to_run','wind_direction','zDep','buoy_loc','lakes','Planet','Model','H0','L0','T','C0','Cg0')
            

    end
    
    pz = plot_height(c,:);
    plot(waveheight_ax,test_speeds(pz ~=0),pz(pz ~=0),'-s','LineWidth',2,'DisplayName',planet_to_run);
    
    legend(waveheight_ax,'show','Location','best')
    drawnow
    

        
end

