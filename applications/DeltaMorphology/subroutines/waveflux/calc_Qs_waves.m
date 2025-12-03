function [Qs_max] = calc_Qs_waves(x,y,wind_mag,wind_angle_deg,rho,rho_s,g)
% Calculates maximum wave-driven littoral transport (Nienhuis+2015) for 
% a shoreline (x,y) undergoing a time series of winds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%   x                 = x coordinates of shoreline [m] [1 x numel(points)] 
%   y                 = y coordinates of shoreline [m] [1 x numel(points)] 
%   wind_mag          = magnitude of wind over a time series [m/s] [1 x numel(timesteps)]
%   wind_angle_deg    = angle of wind (CCW from East) [deg] [1 x numel(timesteps)]
%   rho               = density of liquid [kg/m3] [1 x 1]
%   g                 = gravity of planet [m/s2]  [1 x 1]
% OUTPUTS:
%   Qs_max            = maximum wave-driven littoral transport 
%                       (denominator in eqn. 2, Nienhuis+2015, 
%                       eqn. 6 in suppplement of Nienhuis+2015) [1 x numel(points)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % pre-allocate size  with dummy values
    Qs_max = -999.*size(x);

    make_plot = 0; 
    % Plots:
    % (1) Figure 1: Scatter plot of regional shoreline angles at each point
    % (2) Figure 2: Fetch for all shoreline points for all possible wind directions
    % (3) Figure 3: Fetch for each wind direction in the time series of winds
    % (4) Figure 4: Comparison of non-dimensional growth curve with JONSWAP
    % (5) Figure 5: Plot of shoreline points with labels for point index
    % (6) Figure 6: Wave energy PDFs for each point

    if numel(wind_mag) ~= numel(wind_angle_deg)
        % each value of time has a value of wind_mag and wind_angle_deg
        % not neccessarily all unique values
        % e.g., at time = 1, wind_mag = 1, wind_angle = 0
        %       at time = 2, wind_mag = 1, wind_angle = 0
        %       at time = 3, wind_mag = 2, wind_angle = 10
        error('Calc_Qs_waves: time series of wind speed and angles must be the same size')
    end

    % check that shoreline is a closed vector (beginning = end)
    if abs((x(end) - x(1))) > 1e-6 || abs((y(end) - y(1))) >  1e-6
        disp('Calc_Qs_waves: Shoreline points are not a closed shape, closing the shoreline points by appending start to end')
        x = [x; x(1)];
        y = [y; y(1)];
        added_point = true;
    else
        added_point = false;
    end

    % check if shoreline points are evenly sampled
    all_pts = [x' y']; all_dist = sqrt(sum(diff(all_pts).^2, 2)); tol_dist = 10;
    if max(all_dist) - min(all_dist) > tol_dist
        disp('Points not evenly spaced. Recomputed so points are evenly spaced.')
        % resample along the path with even spacing
        even_spacing = interparc(numel(x), x, y, 'spline'); % https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc
        x = even_spacing(:,1);
        y = even_spacing(:,2);
        x = x(:);
        y = y(:);
    end

    % check that x and y are row vectors
    if ~isrow(x) || ~isrow(y)
        disp('Calc_Qs_waves: Shoreline points given as column vector, converting to row vectors')
        x = x';
        y = y';
    end

     % check that shoreline points are in CCW
    if ispolycw(x,y)
        disp('Calc_Qs_waves: Points given in CW orientation, converting to CCW')
        [x,y] = poly2ccw(x,y);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%
    % Step 1: calculate regional shoreline angle
    %%%%%%%%%

    % shoreline angle is oriented CCW tangential to the chord length
    % connecting points along the window size
    
    window_size_ang = 1;%round(numel(x)/100);  % number of points which angle is being calculated (1% of total number)

    shoreline_angle_rad = calc_regional_shoreline_angle(x, y,window_size_ang); 
    shoreline_angle_deg = rad2deg(shoreline_angle_rad);

    if make_plot
        figure('Name','regional shoreline angle');
        plot(x,y,'-k','LineWidth',2)
        fill(x,y,[0 0.7 0.7],'FaceAlpha',0.2)
        hold on
        scatter(x,y,50,rad2deg(wrapToPi(shoreline_angle_rad)),'filled')
        cb = colorbar();
        clim([-180 180])
        ylabel(cb,'Shoreline Angle [Deg]','FontSize',16, 'Rotation',270)
        axis equal padded
    end
    disp('Regional shoreline angle computed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%
     % Step 2: Calculate fetches for all points in all directions
     %%%%%%%%%
     disp('Fetch calculation begun')
     
     theta_fetch_deg = 0:359; % angles computing fetch over (all possible directions)
     fetch = calc_fetch(x,y,theta_fetch_deg);  

     if make_plot

        figure('Name','Fetch for Each Direction');
        
        for i = 1:numel(theta_fetch_deg)
            clf
            max_fetch = max(fetch(:));
            plot(x,y,'-k','LineWidth',2)
            hold on
            fill(x,y,[0 0.7 0.7],'FaceAlpha',0.2)
            plot_fetch = fetch(i,1:end-1);
            plot_fetch(plot_fetch==0) = NaN;
            scatter(x(1:end-1),y(1:end-1),50,plot_fetch./max_fetch,'filled')
            L = 0.3*(min(range(x),range(y)));
            quiver(mean(x), mean(y), L.*cosd(theta_fetch_deg(i)),  L.*sind(theta_fetch_deg(i)), 1, 'r', 'LineWidth', 2, 'MaxHeadSize', 2); % scale=1
            cb = colorbar();
            ylabel(cb,'Fetch Relative to Max Fetch','FontSize',16, 'Rotation',270)
            clim([0.5 1])
            title(['wind going in direction deg = ', num2str(theta_fetch_deg(i))])
            axis equal padded
            %pause(0.5);
            drawnow;
            if i == 1
                gif('all_fetches.gif')
            else
                gif;
            end   
        end
     end
    disp('Fetch lengths computed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%
    % STEP 3: Get time series of waves at all points
    %%%%%%%%%
    
    % Map each wind direction to its closest directional index
    theta_fetch_deg = theta_fetch_deg(:);
    [~, dir_idx] = min(abs(theta_fetch_deg - wind_angle_deg(:)'), [], 1);  % [1 x numel(wind)]
    fetch_at_pt = fetch(dir_idx, :);                                       %  [numel(wind) x numel(x)]

    if make_plot
        
        figure('Name','Fetch at each Time Step')
        max_fetch = max(fetch(:));

        for i = 1:numel(wind_mag)
            clf
            plot(x,y,'-ok')
            hold on
            fill(x,y,[0 0.7 0.7],'FaceAlpha',0.2)
            plot_fetch = fetch_at_pt(i,1:end-1);
            plot_fetch(plot_fetch==0) = NaN;
            scatter(x(1:end-1),y(1:end-1),100,plot_fetch,'filled')
            L = 0.3*(min(range(x),range(y)));
            quiver(mean(x), mean(y), L.*cosd(wind_angle_deg(i)),  L.*sind(wind_angle_deg(i)), 1, 'r', 'LineWidth', 2, 'MaxHeadSize', 2); % scale=1
            cb = colorbar();
            ylabel(cb,'Fetch [m]','FontSize',16, 'Rotation',270)
            clim([0 max_fetch])
            axis equal padded
            title(['fetches for time t = ',num2str(i)])
            %pause(0.5);
            drawnow;
            if i == 1
                gif('fetch_time_series.gif')
            else
                gif;
            end   
        
        end

    end
    % calculate wave characteristics for a given wind speed, fetch, and gravity
    [wave_height, wave_period, wave_length] = wind2wave(wind_mag, fetch_at_pt, rho,g);  %  [numel(wind) x numel(x)]

    if make_plot

        figure('Name','Wave Height Time Series')

        for  i = 1:numel(x)
            plot(1:numel(wind_mag),wave_height(:,i))
            hold on;
        end

        xlabel('time step')
        ylabel('wave height m')
        title('at each point')
        grid on;
        

    end

    if make_plot
        figure('Name','Wave Height (Rose)')
        for pt = 1:numel(x)
            clf
            plt_H = wave_height(:,pt);
            plt_H(plt_H == 0) = NaN;
            wind_rose(wind_angle_deg,plt_H,'ci',[1 5 10 15])
            title(['Point: ',num2str(pt)])
            drawnow
            if pt == 1
                gif('waveroseforpoints.gif','DelayTime',1)
            else
                gif;
            end
        end
    end

    if make_plot
        % make table comparing estimates for wave height and wave period with JONSWAP
        compare2jonswap(wind_mag,fetch_at_pt,g,wave_height,wave_period)
    end
    disp('Wave time series computed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%
     % STEP 4: Turn time series of waves at each point into Epdf for each point
     %%%%%%%%%
    
     if make_plot
        figure('Name','Point Labels');
        plot(x,y,'-o')
        hold on;
        fill(x,y,[0 0.7 0.7],'FaceAlpha',0.2)
        for pt = 1:numel(x)-1
            text(x(pt),y(pt),num2str(pt),"FontSize",10,'FontWeight','bold')
        end
         hold off
         axis equal padded
     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%
     % STEP 5: Calculate maximum littoral transport of waves
     %%%%%%%%%

     if 1 
        fig = figure('Name','Energy PDFs for Each Point'); 
        t = tiledlayout(fig,1,2);        
        
        % Left panel (scatter)
        axScatter = axes(t);
        axScatter.Layout.Tile = 1;
        
        % Right panel (polar)
        pax = polaraxes(t);
        pax.ThetaZeroLocation = 'right';          % 0 deg to the right 
        pax.ThetaDir = 'counterclockwise';        % make angles increase clockwise
        pax.Layout.Tile = 2;

     end

     for pt = 1:numel(x) 
        
        % convert wind -> wave crest direction
        wave_crests_angle_deg = wrapTo360(wind_angle_deg + 90);         % angle phi0 in Ashton+2006, Figure 1
        
        [wave_crests, Epdf] = make_wave_Epdf(wave_height(:,pt), wave_crests_angle_deg,wave_period(:,pt));
        
       if make_plot

            % Energy PDF
            cla(pax); hold(pax,'on')
            % uncomment next two lines to plot using wind direction
            %polarplot(pax, deg2rad(wrapTo180(wave_crests-90)), Epdf,'LineWidth', 2); % show in direction of where wind is going
            %title(pax, ['Energy PDFs (wind dir) for point: ', num2str(pt)])
            % uncomment next two lines to plot using wave crests
            polarplot(pax, deg2rad(wave_crests), Epdf,'LineWidth', 2); % show in direction of wave crest not wind direction!
            title(pax, ['Energy PDFs (wave crests) for point: ', num2str(pt)])
            rlim(pax, [0 1])
            
            % shoreline
            cla(axScatter); hold(axScatter,'on')
            plot(axScatter, x, y,'-ok')
            scatter(axScatter, x(pt), y(pt), 100, 'r', 'filled') % highlight current pt
            title(axScatter, 'Point Locations')
            xlabel(axScatter, 'x'); ylabel(axScatter, 'y')
            axis(axScatter,'equal')
            axis(axScatter,'padded')
            grid(axScatter,'on')
            
            drawnow;
            
            if pt == 1
                gif('energy_pdfs.gif','DelayTime',1)
            else
                gif;
            end

       end

        Qs_max(pt) = calc_wave_flux_w_PDF(Epdf,wave_height(:,pt),wave_period(:,pt),rho,rho_s,g,wave_crests,shoreline_angle_deg(pt)); % [1 x numel(points)-1] ,wave_height(1:end,1:end-1),wave_period(1:end,1:end-1)
     
        % if pt == numel(x)
        %     Qs_max(pt) = Qs_max(1); % closed shape so start and end are the same
        % end
     end
    
     if added_point
         Qs_max(end) = [];
     end
    disp('Littoral transport computed')
    disp('Wave Qs complete!')
end
