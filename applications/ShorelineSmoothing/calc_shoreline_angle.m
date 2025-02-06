function [wind_dir,wave_front_angle,shoreline_angle,relative_angle] = calc_shoreline_angle(x,y,Wind)
% Initialize arrays for theta and Qs, matching the order of shoreline points (x, y)
relative_angle = NaN(size(x));  % Relative angle of attack

% find theta
for i = 1:numel(x)
    if i < numel(x)
        shoreline_angle(i) = atan2(y(i+1) - y(i), x(i+1) - x(i));
    else % wrap around to the first point to complete the circle
        shoreline_angle(i) = atan2(y(1) - y(i), x(1) - x(i));
    end
    
    shoreline_angle(i) = mod(shoreline_angle(i),2*pi);
    wind_dir = Wind.dir;
    wave_front_angle(i) =  wind_dir + pi / 2;
    relative_angle(i) = mod(wave_front_angle(i) - shoreline_angle(i),2*pi);

    if relative_angle(i) > pi/2 && relative_angle(i) < 3*pi/2
        relative_angle(i) = NaN;
    end
    % Limit angle to the range [0, pi/2]
    % if abs(relative_angle(i)) > pi/2 || abs(relative_angle(i)) < 0
    %     relative_angle(i) = NaN;
    % end
end

end