function [wind_dir, wave_front_angle, shoreline_angle, relative_angle] = calc_shoreline_angle(x, y, Wind,win_size)
    % Initialize arrays for shoreline angle and relative angle
    num_points = numel(x);
    shoreline_angle = NaN(size(x));
    relative_angle = NaN(size(x));

    
    % Calculate shoreline angles using 17 points ahead and behind
    for i = 1:num_points
        % Determine indices for forward and backward points with wrap-around
        idx_forward = mod(i - 1 + win_size, num_points) + 1;
        idx_backward = mod(i - 1 - win_size, num_points) + 1;
        
        % Compute the shoreline angle
        shoreline_angle(i) = atan2(y(idx_forward) - y(idx_backward), x(idx_forward) - x(idx_backward));
        shoreline_angle(i) = mod(shoreline_angle(i), 2 * pi);
    end
    
    % Wind direction and wave front angle
    wind_dir = Wind.dir;
    wave_front_angle = wind_dir + pi / 2;
    
    % Compute relative angles
    for i = 1:num_points
        relative_angle(i) = mod(wave_front_angle - shoreline_angle(i), 2 * pi);
        
    end
end