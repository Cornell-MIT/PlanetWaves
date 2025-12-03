function [x_circle, y_circle] = make_circle(radius,N)

    % Angles from 0 to 360 degrees, step every N degrees
    theta = 0 : N : 360;

    % Convert to radians
    t = deg2rad(theta);

    % Circle coordinates
    x_circle = radius * cos(t);
    y_circle = radius * sin(t);

end