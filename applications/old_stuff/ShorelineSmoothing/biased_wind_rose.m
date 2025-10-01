clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOUSE KEEPING
addpath(fullfile('..','..','planetwaves'))  
addpath(fullfile('..','..','planetwaves/pre_analysis/'))  
%addpath(fullfile('.','waves_w_sloping_bathy'))  

plot_H_grid_chosen = 0;
using_H = 0;

% WAVES FOR 8 DIRECTIONS
wind_direction = 0:45:315;
wind_direction = deg2rad(wind_direction);


window_size_sin = 200;
window_size_ang = 10;

% nice colormaps for plotting
jet_wrap = vertcat(jet,flipud(jet)); % circular colormap
viridis_top = viridis;
viridis_top = viridis_top(1:end/2,:);
plasma_bottom = plasma;
plasma_bottom = plasma_bottom(end/2:end,:);
distinct_cmap = [viridis_top;plasma_bottom]; % distinct colormap to distinguish positive and negative values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(fullfile('.','rednoise_asylake_1')) % not using the files here but still need to reference them
% look here for shoreline, depth, waves (make_deepwater_waves)
fne = {'0deg.mat','45deg.mat','90deg.mat','135deg.mat','180deg.mat','225deg.mat','270deg.mat','315deg.mat'};

shoreline_choice = 3;
% % LAKE SHORELINE WITH SIMPLE BATHYMETRY
if shoreline_choice == 1
    addpath(fullfile('.','rednoise_asylake_1')) % not using the files here but still need to reference them
    load('asylake1.mat','x','y','x_center','y_center')
    % %[x,y] = make_circle(x,y);
    x(end) = [];
    y(end) = [];
    load('lake_depth.mat','Xmesh','Ymesh','zDep')

    figure;
    plot(x,y,'-k')

elseif shoreline_choice == 2
    addpath(fullfile('..','..','\data\Titan\TitanLakes\Bathymetries\bathtub_bathy'))
    load('ol_bathtub_0.46e_neg3.mat','X_cor','Y_cor','x_center','y_center','Xmesh','Ymesh','zDep')
    X_cor(1) = [];
    Y_cor(1) = [];
    x = X_cor;
    y = Y_cor;
    
elseif shoreline_choice == 3
    load('KiskottoStability.mat')
    lon = flip(KiskottoLake.lon);
    lat = flip(KiskottoLake.lat);
    
    lon = rmmissing(lon);
    lat = rmmissing(lat);

    
    wgs84 = geocrs(4326); % EPSG 4326 = WGS84
    
    utmZoneNumber = floor((mean(lon) + 180)/6) + 1;
    isNorth = mean(lat) >= 0;
    
    epsg = 32600 + utmZoneNumber + 100*(~isNorth); 
    utmCrs = projcrs(round(epsg));
    
    [x, y] = projfwd(utmCrs, lat, lon);


    %figure;
    [x,y] = resample_path_even_spacing(x,y,100);
    plot(x,y,'-k')
    hold on
    plot(x(1:window_size_sin),y(1:window_size_sin),'-r')
    title('sinuosity window size')

    figure;
    plot(x,y,'-k')
    hold on
    plot(x(1:window_size_ang),y(1:window_size_ang),'-r')
    title('shoreline angle window size')

    polyin = polyshape(x,y);
    [x_center,y_center] = centroid(polyin);
else
    error('no shoreline, depths provided')
end

% Biased wind rose with a main direction that falls off symmetrically away from main direction as a circular normal
for mu = 1:numel(wind_direction)
    main_direction = wind_direction(mu);
    windPDF = vonMises(main_direction,1.5,wind_direction);
    
    % plot wind pdf
    figure;
    polarplot([wind_direction wind_direction(1)],[windPDF windPDF(1)],'-k','LineWidth',2)
    
    % for each wind direction, calculate wave energy and diffusivity
    for fn = 1:numel(fne)
                
        % calculate sinuosity of the shoreline with window size of 50
        sinu = calc_sinuosity(x,y,window_size_sin);
        
        % calculate relative angle of shoreline with wave crest angle
        Wind.dir = wind_direction(fn);
        Wind.speed = 1;
        % calculate shoreline angle and angle of wavefront relative
        [wind_dir,wave_front_angle,shoreline_angle,relative_angle] = calc_shoreline_angle(x,y,Wind,window_size_ang);
        
         if fn == 1
            % plot sinuosity
            figure;
            scatter(x,y,50,sinu,'filled')
            colorbar
            title('sinuosity')
            % plot shoreline angle (theta)
            figure;
            scatter(x,y,50,wrap_180(rad2deg(shoreline_angle)),'filled')
            colormap(jet_wrap)
            colorbar
            clim([-180 180])
            title('shoreline angle ($\theta$)','Interpreter','latex')
        end
    
        % plot relative angle (psi - theta)
        figure;
        scatter(x,y,50,wrap_180(rad2deg(relative_angle)),"filled")
        colormap(distinct_cmap)
        colorbar
        hold on;
        scatter(x,y,50,'.k')
        quiver(x_center,y_center,1e4*cos(wind_dir),1e4*sin(wind_dir))
        hold off
        title('$\psi$ - $\theta$','Interpreter','latex')
    
        load(fne{fn})

        if using_H == 1
            plot_H = sig_wave.H;
    
            figure;
            h = imagesc(Xmesh(1,:),Ymesh(:,1),plot_H);
            set(h,'AlphaData',~isnan(plot_H))
            axis equal
            xlim([-200 500])
            ylim([-200 500])
            colorbar
            hold on
            plot3(x,y,200.*ones(size(x)),'-r') 
            quiver(x_center,y_center,80*cos(wind_dir),80*sin(wind_dir),'off','MaxHeadSize',5)
            set(gca,'YDir','normal')
            title(sprintf('Waves in %s',fne{fn}))
        end
    
    
        for POI = 1:numel(x)

            if using_H == 1
                % Find nearest Z value
                [H_pt(POI), row, col] = findNearestZ(Xmesh, Ymesh, sig_wave.H, x(POI), y(POI));
                %title(sprintf('wave height is %f', H_pt(POI)));
                T_pt(POI) = sig_wave.T(row,col);
                D_pt(POI) = zDep(row,col);
            end

            H_pt(POI) = 0.9;
            T_pt(POI) = 10.2;
            D_pt(POI) = 80;

            % use deepwater waves
            % [H_pt(POI), i_point] = max(sig_wave.H,[],'all','omitnan');
            % T_pt(POI) = sig_wave.T(i_point);
            % D_pt(POI) = zDep(i_point);

            if plot_H_grid_chosen == 1
                % Create surface plot
                s = imagesc(Xmesh(1,:), Ymesh(:,1), sig_wave.H);
                view(2);
                hold on;
    
                alphaData = ones(size(Xmesh)) * 0.5; 
                alphaData(row, col) = 1; % leep selected grid cell fully visible
                set(s,'AlphaData',alphaData)
    
                 % Plot trajectory and selected point
                plot3(x, y, 100 .* ones(size(x)), '-r');
                plot3(x(POI), y(POI), 100, 'om');
    
                % Highlight selected grid cell
                plot3(Xmesh(row, col), Ymesh(row, col), 100, 'or');
    
                drawnow;
                hold off;
            
                % Save as GIF
                if POI == 1
                    gif('waves.gif');
                else
                    gif;
                end
            end

        end
    
        for ii = 1:numel(x)
            dn_dt(fn,ii) = shoreline_stability(H_pt(ii),T_pt(ii),relative_angle(ii),D_pt(ii));
            wave_energy(fn,ii) = calculate_wave_energy(H_pt(ii),relative_angle(ii));
        end
        weighted_wave_energy(fn,:) = wave_energy(fn,:).*windPDF(fn); 
        weighted_dn_dt(fn,:) = dn_dt(fn,:).*windPDF(fn);
    end


    % figure; 
    % scatter(x,y,50,'.k')
    % hold on
    % scatter(x,y,50,dn_dt(1,:)./max(dn_dt(1,:)),'filled')
    % colorbar

    
    weighted_wave_energy_sum = sum(weighted_wave_energy,1,"omitnan");
    weighted_dn_dt_sum = sum(weighted_dn_dt,1,"omitnan");
    
    weighted_wave_energy_sum = weighted_wave_energy_sum./max(weighted_wave_energy_sum);
    weighted_dn_dt_sum = weighted_dn_dt_sum./max(weighted_dn_dt_sum);
    
    
    
    figure
    scatter(x,y,50,'k')
    hold on
    scatter(x,y,50,weighted_wave_energy_sum,'filled')
    quiver(x_center,y_center,1e4*cos(main_direction),1e4*sin(main_direction))
    colorbar
    clim([0 1])
    title('normalized weighted wave energy with reoccurence')
    
    figure
    scatter(x,y,50,'k')
    hold on
    scatter(x,y,50,weighted_dn_dt_sum,'filled')
    quiver(x_center,y_center,1e4*cos(main_direction),1e4*sin(main_direction))
    colormap(distinct_cmap)
    colorbar
    clim([-1 1])
    title('normalized weighted diffusivity with reoccurence')
    
    [r_wave(mu),s_wave(mu)] = r_squared(weighted_wave_energy_sum,1./sinu,'wave energy','1/sinuoisity');
    [r_diff(mu),s_diff(mu)] = r_squared(weighted_dn_dt_sum,1./sinu,'diffusivity','1/sinuoisity');

    [r_wave_spearman(mu), p_wave_spearman(mu)] = corr(weighted_wave_energy_sum',1./sinu', 'Type', 'Spearman');
    [r_diff_spearman(mu), p_diff_spearman(mu)] = corr(weighted_dn_dt_sum',1./sinu', 'Type', 'Spearman');

    fprintf('For direction %i:\n diff: r = %f p = %f\n wave: r = %f p= %f\n',rad2deg(wind_direction(mu)),r_diff_spearman(mu), p_diff_spearman(mu), r_wave_spearman(mu), p_wave_spearman(mu))


    if mu ~= numel(wind_direction)
        close all
    end

end

% Define the shaded region
theta = linspace(0, 2*pi, 1000); % Full circle angles
r_negative = linspace(-1, 0, 50); % Range of negative r values

% Convert polar to Cartesian for shading
[Theta, R] = meshgrid(theta, r_negative); 
[X, Y] = pol2cart(Theta, R);

% Create figure with polar axes
figure
pax = polaraxes; % Create a polar axes manually
hold on

% Use scatter to simulate shading in the negative r-region
polarscatter(Theta(:), R(:), 10, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.02);

% Plot the main polar plot
polarplot(pax,[wind_direction wind_direction(1)],[s_wave s_wave(1)],'-k','LineWidth',0.5)
polarscatter(pax, wind_direction, s_wave, 50, r_wave, 'filled');
colorbar
% Set radial limits
rlim([-1 1])
hold off
title('slope of wave energy vs 1/sin w color = R^2')

% Create figure with polar axes
figure
pax = polaraxes; % Create a polar axes manually
hold on

% Use scatter to simulate shading in the negative r-region
polarscatter(Theta(:), R(:), 10, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.02);

% Plot the main polar plot
polarplot(pax,[wind_direction wind_direction(1)],[s_diff s_diff(1)],'-k','LineWidth',0.5)
polarscatter(pax, wind_direction, s_diff, 50, r_diff, 'filled');
colorbar
% Set radial limits
rlim([-1 1])
hold off
title('slope of diffusivity vs 1/sin w color = R^2')


p_wave_spearman(p_wave_spearman>0.05) = 1;
p_wave_spearman(p_wave_spearman<= 0.05) = -1;

p_diff_spearman(p_diff_spearman>0.05) = 1;
p_diff_spearman(p_diff_spearman<=0.05) = -1;

wind_direction = wind_direction + pi + pi/2;
figure;
pax2 = polaraxes;
hold on;
polarscatter(pax2,Theta(:), R(:), 10, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.02,'HandleVisibility','off');
polarplot(pax2,[wind_direction wind_direction(1)],[r_wave_spearman r_wave_spearman(1)],'-k','LineWidth',2,'DisplayName','wave energy')
rlim([-1 1])
polarplot(pax2,[wind_direction wind_direction(1)],[r_diff_spearman r_diff_spearman(1)],'--k','LineWidth',2,'DisplayName','diffusivity')
polarscatter(pax2, wind_direction, r_wave_spearman, 50, p_wave_spearman, 'filled','HandleVisibility','off');
polarscatter(pax2, wind_direction, r_diff_spearman, 50, p_diff_spearman, 'filled','HandleVisibility','off');
colormap(flip(coolRedToCoolGreen(10)));
legend('show')


function [x_circle, y_circle] = make_circle(x,y)

    polyin = polyshape(x,y);
    [xc,yc] = centroid(polyin);

    for i = 1:numel(x)
        mydist(i) = sqrt((xc - x(i))^2 + (yc - y(i))^2);
    end

    my_r = max(mydist);

    theta = linspace(0, 2*pi, 1000);
    x_circle = xc + my_r.*cos(theta);
    y_circle = yc + my_r.*sin(theta);

end

function new_angle = wrap_180(angle)
    
    new_angle = mod(angle + 180, 360) - 180;

end

function cmap = coolRedToCoolGreen(n)
    if nargin < 1
        n = 256;
    end

    
    cool_red = [0.8, 0.3, 0.4];  % muted rosy red
    cool_blue = [0.5882, 0.8235,0.5804]; % soft teal blue

    
    r = linspace(cool_red(1), cool_blue(1), n)';
    g = linspace(cool_red(2), cool_blue(2), n)';
    b = linspace(cool_red(3), cool_blue(3), n)';

    cmap = [r g b];
end
