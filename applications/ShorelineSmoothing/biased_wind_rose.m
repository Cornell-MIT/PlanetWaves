clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOUSE KEEPING
addpath(fullfile('..','..','planetwaves'))  
addpath(fullfile('..','..','planetwaves/pre_analysis/'))  
addpath(fullfile('.','waves_w_sloping_bathy'))  

% WAVES FOR 8 DIRECTIONS
wind_direction = 0:45:315;
wind_direction = deg2rad(wind_direction);
fne = {'0deg.mat','45deg.mat','90deg.mat','135deg.mat','180deg.mat','225deg.mat','270deg.mat','315deg.mat'};

% nice colormaps for plotting
jet_wrap = vertcat(jet,flipud(jet)); % circular colormap
viridis_top = viridis;
viridis_top = viridis_top(1:end/2,:);
plasma_bottom = plasma;
plasma_bottom = plasma_bottom(end/2:end,:);
distinct_cmap = [viridis_top;plasma_bottom]; % distinct colormap to distinguish positive and negative values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Biased wind rose with a main direction that falls off symmetrically away from main direction as a circular normal
for mu = 1:numel(wind_direction)
    main_direction = wind_direction(mu);
    windPDF = vonMises(main_direction,1.5,wind_direction);
    
    % plot wind pdf
    figure;
    polarplot([wind_direction wind_direction(1)],[windPDF windPDF(1)],'-k','LineWidth',2)
    
    % for each wind direction, calculate wave energy and diffusivity
    for fn = 1:numel(fne)
        
        % LAKE SHORELINE WITH SIMPLE BATHYMETRY
        load('asylake1.mat')
        load('asylake1_bathtub.mat')
        zDep = zDep.*1000;
        
        zDep(:,1:400) = [];
        zDep(:,180:600) = [];
        zDep(1:440,:) = [];
        zDep(120:560,:) = [];
        
        Xmesh(:,1:400) = [];
        Xmesh(:,180:600) = [];
        Xmesh(1:440,:) = [];
        Xmesh(120:560,:) = [];
        
        Ymesh(:,1:400) = [];
        Ymesh(:,180:600) = [];
        Ymesh(1:440,:) = [];
        Ymesh(120:560,:) = [];
        
        % calculate sinuosity of the shoreline with window size of 100
        sinu = calc_sinuosity(x,y,50);
        % calculate relative angle of shoreline with wave crest angle
        Wind.dir = wind_direction(fn);
        Wind.speed = 1;
        [wind_dir,wave_front_angle,shoreline_angle,relative_angle] = calc_shoreline_angle(x,y,Wind,50);
        
         if fn == 1
            % plot sinuosity
            figure;
            scatter(x,y,50,sinu,'filled')
            colorbar
            title('sinuosity')
            % plot shoreline angle (theta)
            figure;
            scatter(x,y,50,rad2deg(shoreline_angle),'filled')
            colormap(jet_wrap)
            colorbar
            title('shoreline angle (\theta)')
        end
    
        % plot relative angle (psi - theta)
        figure;
        scatter(x,y,50,rad2deg(relative_angle),"filled")
        %colormap(jet_wrap)
        colormap(jet)
        colorbar
        hold on;
        scatter(x,y,50,'.k')
        quiver(x_center,y_center,100*cos(wind_dir),100*sin(wind_dir))
        hold off
        title('\psi - \theta')
    
        % retrieve model parameters
        resizeFactor = 1/10;
        buoy_loc = [200 200];
        grid_resolution = [300 300];
        planet_to_run = 'Titan-OntarioLacus';
        time_to_run = 60;
        Wind.dir = wind_direction(fn);
        [zDep, buoy_loc, grid_resolution, Xmesh, Ymesh, ~, ~] = degrade_depth_mesh(zDep, buoy_loc, grid_resolution, resizeFactor, Xmesh, Ymesh, x, y);
    
    
        load(fne{fn})
        sig_wave = invert_attribute(PeakWave);
    
        figure;
        imagesc(Xmesh(1,:),Ymesh(:,1),sig_wave.H)
        colorbar
        hold on
        plot3(x,y,200.*ones(size(x)),'-r') 
        quiver(x_center,y_center,100*cos(wind_dir),100*sin(wind_dir))
        title(sprintf('Waves in %s',fne{fn}))
    
    
        for POI = 1:numel(x)

            % Find nearest Z value
            [H_pt(POI), row, col] = findNearestZ(Xmesh, Ymesh, sig_wave.H, x(POI), y(POI));
            %title(sprintf('wave height is %f', H_pt(POI)));
            T_pt(POI) = sig_wave.T(row,col);
            D_pt(POI) = zDep(row,col);

            % % Create surface plot
            % s = imagesc(Xmesh(1,:), Ymesh(:,1), sig_wave.H);
            % view(2);
            % hold on;
            % 
            % alphaData = ones(size(Xmesh)) * 0.5; 
            % alphaData(row, col) = 1; % leep selected grid cell fully visible
            % set(s,'AlphaData',alphaData)
            % 
            %  % Plot trajectory and selected point
            % plot3(x, y, 100 .* ones(size(x)), '-r');
            % plot3(x(POI), y(POI), 100, 'om');
            % 
            % % Highlight selected grid cell
            % plot3(Xmesh(row, col), Ymesh(row, col), 100, 'or');
            % 
            % drawnow;
            % hold off;
        
            % % Save as GIF
            % if POI == 1
            %     gif('waves.gif');
            % else
            %     gif;
            % end
            % 
        end
    
    
    
        for ii = 1:numel(x)
            dn_dt(fn,ii) = shoreline_stability(H_pt(ii),T_pt(ii),relative_angle(ii),D_pt(ii));
            wave_energy(fn,ii) = H_pt(ii).^2.*cos(relative_angle(ii));
            if cos(relative_angle(ii)) < 0
                wave_energy(fn,ii) = NaN;
            end
        end
        weighted_wave_energy(fn,:) = wave_energy(fn,:).*windPDF(fn); 
        weighted_dn_dt(fn,:) = dn_dt(fn,:).*windPDF(fn);
    end
    
    
    weighted_wave_energy_sum = sum(weighted_wave_energy,1,"omitnan");
    weighted_dn_dt_sum = sum(weighted_dn_dt,1,"omitnan");
    
    weighted_wave_energy_sum = weighted_wave_energy_sum./max(weighted_wave_energy_sum);
    weighted_dn_dt_sum = weighted_dn_dt_sum./max(weighted_dn_dt_sum);
    
    
    
    
    
    figure
    scatter(x,y,50,'k')
    hold on
    scatter(x,y,50,weighted_wave_energy_sum,'filled')
    quiver(x_center,y_center,100*cos(main_direction),100*sin(main_direction))
    colorbar
    clim([0 1])
    title('normalized weighted wave energy with reoccurence')
    
    figure
    scatter(x,y,50,'k')
    hold on
    scatter(x,y,50,weighted_dn_dt_sum,'filled')
    quiver(x_center,y_center,100*cos(main_direction),100*sin(main_direction))
    colormap(distinct_cmap)
    colorbar
    clim([-1 1])
    title('normalized weighted diffusivity with reoccurence')
    
    [r_wave(mu),s_wave(mu)] = r_squared(weighted_wave_energy_sum,1./sinu,'wave energy','1/sinuoisity');
    [r_diff(mu),s_diff(mu)] = r_squared(weighted_dn_dt_sum,1./sinu,'diffusivity','1/sinuoisity');
    close all

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
polarscatter(Theta(:), R(:), 10, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.5);

% Plot the main polar plot
polarplot(pax, [wind_direction wind_direction(1)], ...
              [r_wave./sign(s_wave) r_wave(1)/sign(s_wave(1))], ...
              '-k', 'LineWidth', 2);

% Set radial limits
rlim([-1 1])

hold off
title('R^2/sign(slope) wave energy vs 1/sin')


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
polarscatter(Theta(:), R(:), 10, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.5);

% Plot the main polar plot
polarplot(pax, [wind_direction wind_direction(1)], ...
              [r_diff./sign(s_diff) r_diff(1)/sign(s_diff(1))], ...
              '-k', 'LineWidth', 2);

% Set radial limits
rlim([-1 1])

hold off
title('R^2/sign(slope) diffusivity vs 1/sin')