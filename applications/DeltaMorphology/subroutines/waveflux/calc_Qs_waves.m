function [Qs_total,Qs_time] = calc_Qs_waves(x,y,wind_mag,wind_angle_deg,rho,g)
% angle in degrees (degrees CCW from east)
    
    make_plot = 0;

   
    if numel(wind_mag) ~= numel(wind_angle_deg)
        % each value of time has a value of wind_mag and wind_angle_deg
        % not neccessarily all unique values
        % e.g., at time = 1, wind_mag = 1, wind_angle = 0
        %       at time = 2, wind_mag = 1, wind_angle = 0
        %       at time = 3, wind_mag = 2, wind_angle = 10
        error('Calc_Qs_waves: time series of wind speed and angles must be the same size')
    end

    if ispolycw(x,y)
        disp('Points given in CW orientation, resorting to CCW')
        [x,y] = poly2ccw(x,y);
    end


    Qs_total = NaN;
    Qs_time = NaN;

    %%%%%%%%%
    % Step 1: calculate shoreline angle
    %%%%%%%%%

    window_size_ang = round(numel(x)/5); 
    shoreline_angle_rad = calc_regional_shoreline_angle(x, y,window_size_ang);
    sang = rad2deg(shoreline_angle_rad);

    if make_plot
        figure('Name','shoreline angle');
        plot(x,y,'-k','LineWidth',2)
        hold on
        scatter(x,y,50,rad2deg(wrapToPi(shoreline_angle_rad)),'filled')
        colorbar
    end

     thetas = -90:90;
     Qsmax = zeros(size(sang));
 
     %%%%%%%%%
     % Step 2: Calculate fetches for all points in all directions
     %%%%%%%%%

     theta_fetch_deg = 0:359;
     fetch_km = calc_fetch(x,y,theta_fetch_deg);
     fetch = fetch_km.*1000;

     
     if make_plot
        figure('Name','Max Fetch');
        for i = 1:numel(theta_fetch_deg)
            clf
            max_pos = sqrt((max(x)*1000)^2 + (max(y)*1000)^2);
            plot(x,y,'-k','LineWidth',2)
            hold on
            plot_fetch = fetch(i,1:end-1);
            plot_fetch(plot_fetch==0) = NaN;
            scatter(x(1:end-1),y(1:end-1),50,plot_fetch./max_pos,'filled')
            text(mean(x),mean(y),num2str(round(plot_fetch(~isnan(plot_fetch)))./max_pos))
            quiver(mean(x), mean(y), mean(x).*cosd(theta_fetch_deg(i)),  mean(x).*sind(theta_fetch_deg(i)), 1, 'r', 'LineWidth', 2, 'MaxHeadSize', 2); % scale=1
            colorbar
            clim([0.5 1])
            axis([min(x)-1000 max(x)+1000 min(y)-1000 max(y)+1000]);
            title(['deg = ', num2str(theta_fetch_deg(i))])
            pause(0.5);
            % drawnow;
            % if i == 1
            %     gif('all_winds.gif')
            % else
            %     gif;
            % end
            
        end
     end

    %%%%%%%%%
    % STEP 3: Get time series of waves at all points
    %%%%%%%%%

    wave_height = NaN(numel(x), numel(wind_mag));
    wave_period = NaN(numel(x), numel(wind_mag));
    wave_length = NaN(numel(x), numel(wind_mag));
    
    
    % Map each wind direction to its closest directional index
    theta_fetch_deg = theta_fetch_deg(:);
    [~, dir_idx] = min(abs(theta_fetch_deg - wind_angle_deg(:)'), [], 1);  % [1 x numel(wind)]
    fetch_at_pt = fetch(dir_idx, :);   %  [numel(wind) x numel(x)]

    if 1
        
        figure('Name','Fetch at each Time Step')
        
        for i = 1:numel(wind_mag)
            clf
            plot(x,y,'-ok')
            hold on
            plot_fetch = fetch_at_pt(i,1:end-1);
            plot_fetch(plot_fetch==0) = NaN;
            scatter(x(1:end-1),y(1:end-1),100,plot_fetch,'filled')
            quiver(mean(x), mean(y), mean(x).*cosd(wind_angle_deg(i)),  mean(x).*sind(wind_angle_deg(i)), 1, 'r', 'LineWidth', 2, 'MaxHeadSize', 2); % scale=1
            colorbar
            axis([min(x)-1000 max(x)+1000 min(y)-1000 max(y)+1000]);
            title(['t = ',num2str(i)])
        
        end

    end
    % calculate wave characteristics for a given wind speed, fetch, and gravity
    [wave_height, wave_period, wave_length] = wind2wave(wind_mag, fetch_at_pt, g);  %  [numel(wind) x numel(x)]

    if make_plot
        % make table comparing estimates for wave height and wave period with JONSWAP
        compare2jonswap(wind_mag,fetch_at_pt,g,wave_height,wave_period)
    end

     %%%%%%%%%
     % STEP 4: Turn time series of waves at each point into Epdf for each point
     %%%%%%%%%
    
     if 1
        figure('Name','Point Labels');
        plot(x,y,'-o')
        hold on;
        for pt = 1:numel(x)-1
            text(x(pt),y(pt),num2str(pt))
        end
         hold off
         axis equal
     end

     for pt = 1:numel(x)-1
        [theta_bins, Epdf(pt,:)] = make_wave_Epdf(wave_height(:,pt), wind_angle_deg,rho,g);
        title(['For point: ',num2str(pt)])
        hold off;
    end
     
    
end

 


