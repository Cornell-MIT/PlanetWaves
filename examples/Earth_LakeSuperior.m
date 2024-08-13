clc
clear
close all

% LAKE SUPERIOR WITH REAL BUT DEGRADED BATHYMETRY

warning('need to streamline process of moving outputs of find_fetch.py to matlab scripts to run model')

addpath(fullfile('..','planetwaves'))  
addpath(fullfile('..','planetwaves','pre_analysis'))
addpath(fullfile('..','data','Earth','GreatLakes','LakeSuperior'))
% from find_fetch.py 
load('..\data\Earth\GreatLakes\LakeSuperior\BathyData\LakeSuperior_cleaned.mat')
zDep = -squeeze(LS);

% MODEL INPUTS
planet_to_run = 'Earth';
buoy_loc = [6618 1729];                                                    % grid location [x,y]
% from find_fetch.py
grid_resolution = [4542.948547909539 92.66280063299297];                   % pixel width and pixel height [m]
test_speeds = 0:20;                                                        % wind speed
time_to_run = 60*10;                                                       % time to run model
wind_direction = 0;                                                        % wind direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL MODEL
% degrade depth profile so model doesnt take as long to run
[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.02);
zDep = max(max(zDep)).*ones(12,12);
% populate model classes
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
% adjust anenometer height to buoy
Model.z_data = 3.6;

% update grid resolution
Model.gridX = 20*1000;                                              
Model.gridY = 20*1000;                                               

Model.long = 6;                                                  
Model.lat = 6; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN MODEL

% Preallocate cell arrays to store results
myHsig = cell(1, numel(test_speeds));
htgrid = cell(1, numel(test_speeds));
E_spec = cell(1, numel(test_speeds));
Cg = cell(1,numel(test_speeds));

for i = 1:numel(test_speeds)

    Wind.speed = test_speeds(i);

    [myHsig{i}, htgrid{i}, wn_e_spectrum, ~ , ~ , ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
    if ~isempty(wn_e_spectrum{end})
        energy{i} = squeeze(sum(wn_e_spectrum{end}.E(Model.long,Model.lat,:,:),4));
        wn{i} = squeeze(sum(wn_e_spectrum{end}.k(Model.long,Model.lat,:,:),4));
        cg{i} = squeeze(sum(wn_e_spectrum{end}.cg(Model.long,Model.lat,:,:),4));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PLOT RESULTS
save('lakesuperior_run.mat','myHsig','energy','wn','cg')

figure
for i = 1:numel(test_speeds)
    plot(myHsig{i}, '-', 'LineWidth', 3, 'DisplayName', num2str(test_speeds(i)))
    hold on
end
grid on;
legend('show', 'Location', 'northwest','interpreter','latex');
title(['Waves on',' ',Planet.name],'interpreter','latex');
xlabel('model time step [$\Delta$ t]','interpreter','latex')
ylabel('significant wave height [m]','interpreter','latex')


for k = 1:numel(test_speeds)
    buoy_waves(k) = htgrid{k}{end}(Model.long,Model.lat);
end

figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1])
for speed = 1:numel(test_speeds)

    subplot(1,3,3)
    plot(test_speeds,buoy_waves,'-sb','LineWidth',1,'MarkerFaceColor','b')
    hold on
    plot(test_speeds(speed),htgrid{speed}{end}(Model.long,Model.lat),'-sr','LineWidth',1,'MarkerFaceColor','r')
    hold off;
    xlabel('u [m/s]','interpreter','latex')
    ylabel('Hsig [m]','interpreter','latex')
    grid on;

    plot_grid = htgrid{speed}{end}';
    plot_grid(isnan(plot_grid)) = 0;
    plot_alpha_data = ones(size(plot_grid));
    plot_alpha_data(plot_grid==0) = 0;

    ax1 = subplot(1,3,[1,2]);
    h1 = imagesc(plot_grid);
    colormap linspecer
    xlabel('longitude [km]')
    ylabel('latitude [km]')
    title(sprintf('u = %i m/s',test_speeds(speed)))
    c1 = colorbar;
    c1.Label.String = 'Hsig [m]';
    clim([0 6])
    set(h1, 'AlphaData', plot_alpha_data);
    hold on;

    contour(Model.bathy_map,linspace(min(min(Model.bathy_map)),max(max(Model.bathy_map)),10),'-k','LineWidth',2)

    grid on
    new_xtick = get(gca, 'XTick')*(Model.gridX)/1000;
    new_ytick = get(gca, 'YTick')*(Model.gridY)/1000;
    set(gca, 'XTick',  get(gca, 'XTick'), 'XTickLabel', arrayfun(@(x) sprintf('%d', x), new_xtick, 'UniformOutput', false));
    set(gca, 'YTick',  get(gca, 'YTick'), 'YTickLabel', arrayfun(@(y) sprintf('%d', y), new_ytick, 'UniformOutput', false));


    [wx,wy] = pol2cart(Wind.dir,1);
    plot(Model.long,Model.lat,'pentagram','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20)

    for i = 1:Model.LonDim
        for j = 1:Model.LatDim
            if Model.bathy_map(j,i) <=0
                quiver(i, j, (speed/5)*wx, (speed/5)*wy, 'k', 'MaxHeadSize', 1);
            end
        end
    end
    set(ax1,'Ydir','reverse')
    %set(ax1,'Xdir','reverse')
    drawnow
    if speed == 1
        gif('LakeSuperior_Uniform.gif','DelayTime',1,'overwrite',true)
    else
        gif
    end
end



figure('units','normalized','outerposition',[0 0 1 1])
plot(test_speeds,buoy_waves,'-sb','LineWidth',1,'MarkerFaceColor','b')
ylim([0 6])
xlabel('u [m/s]','interpreter','latex')
ylabel('Hsig [m]','interpreter','latex')
grid on;

%%%



f_vec=(log(Model.max_freq)-log(Model.min_freq))/(Model.Fdim-1);             
f = exp(log(Model.min_freq)+(0:Model.max_freq-1)*f_vec);  
dang =  360/Model.Dirdim; 

figure('units', 'normalized', 'outerposition', [0 0 1 1]);


color1 = [229, 211, 82] / 255;
color2 = [172, 57, 49] / 255;

% Number of colors in the colormap
num_colors = numel(test_speeds);

% Generate the colormap
cmap = [linspace(color1(1), color2(1), num_colors)', ...
                   linspace(color1(2), color2(2), num_colors)', ...
                   linspace(color1(3), color2(3), num_colors)'];


set(gca,'ColorOrder',colormap(cmap))
hold on;


for i = 1:length(test_speeds)
    
    k_wave = wn{i};
    e_wave = energy{i}*dang;
    plot(k_wave,e_wave,'-','DisplayName',num2str(test_speeds(i)),'LineWidth',5)
    hold on;

    [~,max_ind] = max(e_wave);
    inds = k_wave >= k_wave(max_ind);

    X = k_wave(inds); 
    logX = log(X);
    Y = e_wave(inds); 
    logY = log(Y);
  
    validInds = ~isinf(logY);
    logX = logX(validInds);
    X = X(validInds);
    Y = Y(validInds);
    logY = logY(validInds);

    coeffs = polyfit(logX, logY, 1);
    n = abs(coeffs(1));
    fittedLogY = polyval(coeffs, logX);

    residuals = logY - fittedLogY;
    stdResiduals = std(residuals);
    Sxx = sum((logX - mean(logX)).^2);
    n_std = stdResiduals / sqrt(Sxx);

    n_loop(i) = n;
    n_std_loop(i) = n_std(1);

end



meanY = X.^(-n)*exp(coeffs(2));
plot(X, meanY, '--k', 'LineWidth', 3, 'DisplayName', 'fit')

legend('show')
grid on;
xlabel('k [m-1]', 'FontSize', 25, 'interpreter', 'latex')
ylabel('E(k)', 'FontSize', 25, 'interpreter', 'latex')
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
ylim([1e-6 1e5])
set(gca, 'FontSize', 16)
set(gca, 'FontWeight', 'bold')
box on;


title(sprintf('High frequency tail: f^n where n is: %f ± %f\n', mean(n_loop), mean(n_std_loop)))