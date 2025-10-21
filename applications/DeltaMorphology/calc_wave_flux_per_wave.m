function [Qsmax_total, Qsmax_ts] = calc_wave_flux_per_wave(wind_angle, wind_mag, x,y,sang)
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
    window_size_ang = 50; % size of window to define the regional shoreline orientation over
    shoreline_angle_rad = calc_regional_shoreline_angle(x, y,window_size_ang);
    sang = rad2deg(shoreline_angle_rad);

    thetas = -90:90;                        % Delta orientations relative to shoreline
    T = numel(wind_angle);                  % Number of time steps
    N = numel(sang);                        % Number of shoreline angles

    Qsmax_total = zeros(1, N);              % Output: total Qsmax per shoreline
    Qsmax_ts = zeros(N, T);                 % Output: per-time-step Qsmax per shoreline

    % Loop over shoreline orientations
    for c = 1:N
        regional_shoreline = sang(c);
        Qsnet_total = zeros(size(thetas));  % Accumulate over time

        for t = 1:T
            angle = wind_angle(t);

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

        % Compute total Qsmax after all time steps
        Qsmaxright_total = max(Qsnet_total(thetas <= 0));
        Qsminleft_total = min(Qsnet_total(thetas >= 0));
        Qsmax_total(c) = Qsmaxright_total - Qsminleft_total;
    end
end
