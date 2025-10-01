function shoreline_angle = calc_shoreline_ang(x,y,win_size)
% angle output as radians

% Calculate shoreline angles using window size points ahead and behind
    for i = 1:numel(x)
        idx_forward = mod(i - 1 + win_size, numel(x)) + 1;
        idx_backward = mod(i - 1 - win_size, numel(x)) + 1;

        shoreline_angle(i) = atan2(y(idx_forward) - y(idx_backward), x(idx_forward) - x(idx_backward));
        shoreline_angle(i) = mod(shoreline_angle(i), 2 * pi);
    end
end

