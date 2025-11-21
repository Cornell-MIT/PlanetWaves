function [x,y] = make_synthetic_shoreline(N)


    rng(42); % random seed
    angles = linspace(0, 2*pi, N);


    radii = 5 + 0.2.*randn(1, N) + 5.*sin(3*angles);
    radii = radii.*1000;
    x = radii .* cos(angles);
    y = radii .* sin(angles);
    x = [x x(1)];
    y = [y y(1)];

end