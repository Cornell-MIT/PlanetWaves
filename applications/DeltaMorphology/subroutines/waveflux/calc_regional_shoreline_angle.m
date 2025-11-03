function shoreline_angle = calc_regional_shoreline_angle(x, y,window_size)
    
    
    make_plot = 0;
    num_points = numel(x); shoreline_angle = NaN(size(x)); 

    % force CCW so normal points in to lake and not out from lake
    [x,y] = make_CCW(x,y);
    
    % Calculate shoreline angles using window size points ahead and behind
    for i = 1:num_points
        idx_forward = mod(i - 1 + window_size, num_points) + 1;
        idx_backward = mod(i - 1 - window_size, num_points) + 1;

        shoreline_angle(i) = atan2(y(idx_forward) - y(idx_backward), x(idx_forward) - x(idx_backward));
        shoreline_angle(i) = mod(shoreline_angle(i), 2 * pi);
    end
    
    % Compute normal direction (pointing to liquid side for CCW shoreline)
    shoreline_angle = mod(shoreline_angle + pi/2, 2 * pi); % defining angle of the normal vector pointing to liquid
   
    if make_plot
           figure; hold on; axis equal;
           plot(x, y, 'k-', 'LineWidth', 1.5);
           hold on
           scale_arrow = 1;
           u = cos(shoreline_angle);
           v = sin(shoreline_angle);
           quiver(x, y, scale_arrow.*u, scale_arrow.*v, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);

    end
   
end