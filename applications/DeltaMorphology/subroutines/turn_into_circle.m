function [xc,yc] = turn_into_circle(x,y)


    x = x(:);
    y = y(:);

    % Center of the circle
    xc0 = mean(x);
    yc0 = mean(y);

    % Radius based on data range
    r = max(range(x), range(y)) / 2;

    theta = deg2rad(0:10:360);

    % Circle coordinates
    xc = xc0 + r * cos(theta);
    yc = yc0 + r * sin(theta);
end

