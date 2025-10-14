function shoreline_angle = calc_regional_shoreline_angle(x, y,window_size)
    
    num_points = numel(x); shoreline_angle = NaN(size(x)); 
    
    % Calculate shoreline angles using window size points ahead and behind
    for i = 1:num_points
        idx_forward = mod(i - 1 + window_size, num_points) + 1;
        idx_backward = mod(i - 1 - window_size, num_points) + 1;

        shoreline_angle(i) = atan2(y(idx_forward) - y(idx_backward), x(idx_forward) - x(idx_backward));
        shoreline_angle(i) = mod(shoreline_angle(i), 2 * pi);
    end
    
    shoreline_angle = mod(shoreline_angle + pi/2, 2 * pi); % defining angle of the normal vector pointing to liquid
    
   
end