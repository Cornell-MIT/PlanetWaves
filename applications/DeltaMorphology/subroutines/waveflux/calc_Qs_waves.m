function [Qs_total,Qs_time] = calc_Qs_waves(x,y,wind_mag,wind_angle_deg,rho,g)
% angle in degrees (degrees CCW from east)
    
    make_plot = 0;
        
    
    % each value of time has a value of wind_mag and wind_angle_deg
    % not neccessarily all unique values
    % e.g., at time = 1, wind_mag = 1, wind_angle = 0
    %       at time = 2, wind_mag = 1, wind_angle = 0
    %       at time = 3, wind_mag = 2, wind_angle = 10

    if numel(wind_mag) ~= numel(wind_angle_deg)
        error('Calc_Qs_waves: time series of wind speed and angles must be the same size')
    end

    if ispolycw(x,y)
        disp('Points given in CW orientation, resorting to CCW')
        [x,y] = poly2ccw(x,y);
    end

    time = 1:numel(wind_mag); 

    Qs_total = NaN;
    Qs_time = NaN;

    

    % Step 1: calculate shoreline angle
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
 
     % Step 2: Calculate fetches for all points in all directions
     theta_fetch_deg = 0:359;
     fetch_km = calc_fetch(x,y,theta_fetch_deg);
     fetch = fetch_km.*1000;

     
     if make_plot
        figure('Name','Max Fetch');
        plot(x,y,'-k','LineWidth',2)
        hold on
        scatter(x,y,50,max(fetch),'filled')
        colorbar
     end

     % STEP 3: Get time series of waves at all points

     % loop over time
     for t = 1:numel(time)
         u = wind_mag(t);
         u_ang_deg = wind_angle_deg(t);
         % loop over all shoreline points
         for pt = 1:numel(x)

             fetch_at_pt = fetch(:,pt); % fetch at a point organized 0 - 359 deg
             i_ang = find(theta_fetch_deg == u_ang_deg);
             % get time series of wave properties at point
             [wave_height(pt,t),wave_period(pt,t),wave_length(pt,t)] = wind2wave(u,fetch_at_pt(i_ang),g);

         end

         if make_plot
             make_lake_wind_plot(x,y,wind_angle_deg(t));
             scatter(x,y,50,wave_height(:,t),'filled')
             colorbar
             title(['wave heights at time = ',num2str(t)])
             hold off
         end

     end

     
     % STEP 4: Turn time series of waves at each point into Epdf for each point
    
    if 1
        figure('Name','Point Labels');
        plot(x,y,'-o')
        hold on;
        for pt = 1:numel(x)-1
            text(x(pt),y(pt),num2str(pt))
        end
        hold off
    end

     for pt = 1:numel(x)-1
        [theta_bins, Epdf(pt,:)] = make_wave_Epdf(wave_height(pt,:), wind_angle_deg,rho,g);
        title(['For point: ',num2str(pt)])
        hold off;


    end

    for t = 1:numel(wind_mag)
        ang = find(theta_fetch_deg == wind_angle_deg(t));
        plot_fetch = fetch(ang,:);
        plot_fetch(plot_fetch==0) = NaN;
        figure('Name','Fetch');
        plot(x,y,'-k','LineWidth',2)
        hold on
        scatter(x,y,50,plot_fetch,'filled')
        colorbar
        title(['For time: ',num2str(t)])
        axis equal;
        hold off;
        
    end
     
    
end

 


