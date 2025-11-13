function [x_circle, y_circle] = make_circle(x,y)

% Step 1: Calculate the centroid
centroid_x = mean(x);
centroid_y = mean(y);

% Step 2: Calculate the maximum distance from the centroid to each point
distances = sqrt((x - centroid_x).^2 + (y - centroid_y).^2);
radius = min(distances)/5;

% Step 3: Create a circle with this centroid and radius
theta = linspace(0, 2*pi, 100);  % 100 points to approximate the circle
x_circle = centroid_x + radius * cos(theta);
y_circle = centroid_y + radius * sin(theta);

x_circle(end) = [];
y_circle(end) = [];

end