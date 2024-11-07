function theta = calc_shoreline_angle(x,y,Wind)
% Initialize arrays for theta and Qs, matching the order of shoreline points (x, y)
theta = NaN(size(x));  % Relative angle of attack

% find theta
for i = 1:numel(x)
    if i < numel(x)
        shoreline_angle(i) = atan2(y(i+1) - y(i), x(i+1) - x(i));
    else % wrap around to the first point to complete the circle
        shoreline_angle(i) = atan2(y(1) - y(i), x(1) - x(i));
    end
    
    shoreline_angle(i) = wrapToPi(shoreline_angle(i));
    % Calculate the relative wave angle
    wave_front_angle = wrapToPi(Wind.dir + pi / 2);
    theta(i) = wrapToPi(wave_front_angle - shoreline_angle(i));

    % % Limit theta to the range [0, pi/2]
    if abs(theta(i)) > pi/2 || abs(theta(i)) < 0
        theta(i) = NaN;
    end
end

end