function [x_circle, y_circle] = make_circle(x,y,N)

    polyin = polyshape(x,y);
    [xc,yc] = centroid(polyin);

    for i = 1:numel(x)
        mydist(i) = sqrt((xc - x(i))^2 + (yc - y(i))^2);
    end

    my_r = max(mydist);

    theta = linspace(0, 2*pi, N);
    x_circle = xc + my_r.*cos(theta);
    y_circle = yc + my_r.*sin(theta);

end