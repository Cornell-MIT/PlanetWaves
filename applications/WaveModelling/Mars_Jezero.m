clc
clear
close all

% PLOT THE WAVES IN JEZERO CRATER LAKE 


addpath(genpath(fullfile('..','..','planetwaves')))
addpath(fullfile('..','..','data','Mars'))

% Model conditions at Mars
planet_to_run = 'Mars-low';
time_to_run = 60*10;
test_speeds = 2:0.1:3;%[1.7 2:10];
wind_direction = pi;

% Jezero crater DTM
fn = 'M20_JezeroCrater_CTXDEM_20m.tif';
lake_level = 10;
A = read(Tiff(fn,'r'));
% clean up Tiff
A(A<-1e10) = NaN;
A(A>-2400) = 0;
A(A>0) = 0;
% isolate crater lake
A(:,1:500) = [];
A(:,2900:end) = [];
A(1:1500,:) = [];
A(2500:end,:) = [];
grid_resolution = [20 20]; % 20 m/pix
buoy_loc = [580 836]; % location of Jezero delta

% degrade depth resolution to make it not huge data-wise
[A,buoy_loc,grid_resolution] = degrade_depth_resolution(A,buoy_loc,grid_resolution,0.01);

% setting water level within the crater
A(A>0) = 0;
min_original = min(min(A)); 
max_original = max(max(A));
min_new = -lake_level;
max_new = 0;
A = (A - min_original) / (max_original - min_original) * (max_new - min_new) + min_new;
A = -A;
A(A<1) = 0;


[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,A,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);  

make_input_map(Planet,Model,Wind)

figure('Name',['Waves on ',planet_to_run]);
hold on

for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);
    Model = calc_cutoff_freq(Planet,Model,Wind);

    [myHsig, htgrid, ~, ~ , ~ , ~, PeakWaves] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
     % if ~isempty(wn_e_spectrum{end})
     %    energy{i} = squeeze(sum(wn_e_spectrum{end}.E(Model.long,Model.lat,:,:),4));
     %    wn{i} = squeeze(sum(wn_e_spectrum{end}.k(Model.long,Model.lat,:,:),4));
     %    cg{i} = squeeze(sum(wn_e_spectrum{end}.cg(Model.long,Model.lat,:,:),4));
     % end

    T_p = PeakWaves.T_weighted;
    if ~isempty(htgrid{end})
        H_low(i) = htgrid{end}(Model.long,Model.lat);
        T_low(i) = T_p(Model.long,Model.lat);
    else
        H_low(i) = NaN;
        T_low(i) = NaN;
    end
    yyaxis right
    plot(test_speeds(i),T_low(i),'-sr','MarkerFaceColor','r','LineWidth',5,'MarkerSize',15,'DisplayName', 'low')
    ylabel('Wave Period, T [s]')
    yyaxis left
    plot(test_speeds(i),H_low(i),'-sk','MarkerFaceColor','k','LineWidth',5,'MarkerSize',15,'DisplayName', 'low')
    ylabel('Signifigant Wave Height, H_s [m]')
    drawnow;



end


% save('MarsJezero_low.mat','myHsig','htgrid','T_p')
% make_plots(Planet,Model,Wind,test_speeds,myHsig,htgrid,energy,wn)