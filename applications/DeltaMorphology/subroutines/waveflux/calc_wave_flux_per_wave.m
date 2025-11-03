function [Qsmax_total, Qsmax_ts] = calc_wave_flux_per_wave(wind_angle, wind_mag, x,y)
    % Calculates Qsmax over a time series of wave angle/mag for multiple shoreline orientations
    % 
    % INPUTS:
    %   wind_angle: [1 x T] wave directions (degrees CCW from East, coming from)
    %   wind_mag:   [1 x T] wave magnitudes (e.g., energy or height)
    %   x,y:        [1 x N] shoreline points
    %
    % OUTPUTS:
    %   Qsmax_total: [1 x N] total Qsmax per shoreline over time
    %   Qsmax_ts:    [N x T] time series of Qsmax per shoreline

    % Calculate the regional shoreline orientaiton at each point
    window_size_ang = 10; % size of window to define the regional shoreline orientation over
    shoreline_angle_rad = calc_regional_shoreline_angle(x, y,window_size_ang);
    sang = rad2deg(shoreline_angle_rad);

    thetas = -90:90;                        % Delta orientations relative to shoreline
    T = numel(wind_angle);                  % Number of time steps
    N = numel(sang);                        % Number of shoreline angles

    Qsmax_total = zeros(1, N);              % Output: total Qsmax per shoreline
    Qsmax_ts = NaN(N, T);                 % Output: per-time-step Qsmax per shoreline
    
    figure;
    Qsnet_total = zeros(size(thetas));  % Accumulate over time is this keep resetting to zero?

    for t = 1:T     % Loop over time
        angle = wind_angle(t);
        for c = 1:N % loop over shoreline points
            regional_shoreline = sang(c);
            

            fetch = calc_fetch(x,y,c,mod(angle + 180, 360)); 
            [wave_height,wave_period,~] = wind2wave(wind_mag(t),fetch);
            phi0 = wrapTo180(angle - regional_shoreline);

            if phi0 < -90 || phi0 > 90
                Qsmax_ts(c, t) = 0;  % Wave from land — no transport
                continue;
            end

            Qsnet = zeros(size(thetas));
            for i = 1:numel(thetas)
                theta = thetas(i);
                LST = CERC(wave_height, wave_period, wrapTo180(phi0 - theta));
                warning('calc_wave_flux_per_wave: need to replace CERC with Diegaard')
                if isnan(LST)
                    LST = 0; 
                end
                Qsnet(i) = LST;
                Qsnet_total(i) = Qsnet_total(i) + LST;  % accumulate for total
            end

            % Compute per-time-step Qsmax
            Qsmaxright = max(Qsnet(thetas <= 0));
            Qsminleft = min(Qsnet(thetas >= 0));
            Qsmax_ts(c, t) = Qsmaxright - Qsminleft;
          
        end

            clf
            hold on;
            axis equal;
            
            % Plot shoreline polygon
            plot(x, y, 'k-', 'LineWidth', 2);
            scatter(x, y, 50,Qsmax_ts(:,t), 'filled' );
                       
            % Plot wind direction arrows
            % Convert wind angle (degrees CCW from East) to dx, dy for quiver
            theta_rad = deg2rad(wind_angle(1));  % constant wind
            dx = cos(theta_rad);
            dy = sin(theta_rad);
            
            % Scale arrows
            arrow_scale = 10;
            
            quiver(mean(x), mean(y), dx*arrow_scale, dy*arrow_scale, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
            text(mean(x)+dx*arrow_scale, mean(y)+dy*arrow_scale, 'Wind', 'Color','r','FontSize',12);
            
            xlabel('X');
            ylabel('Y');
            title('Shoreline and Wind Direction');
            grid on;
            hold off;
            drawnow

    end
            % Compute total Qsmax after all time steps
        Qsmaxright_total = max(Qsnet_total(thetas <= 0));
        Qsminleft_total = min(Qsnet_total(thetas >= 0));
        Qsmax_total(c) = Qsmaxright_total - Qsminleft_total;
end
