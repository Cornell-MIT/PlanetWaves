function [Xmesh, Ymesh, zDep] = make_bathtub_lake_w_azimuth_asymmetry(bath_slope, shoreline)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective:
%   This function will create the mesh for a bathtub model of lake bathymetry.
%   Each point on the shoreline has a unique slope defined in bath_slope.
%   The contour3 function requires the bathymetry to be in mesh form to work.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   bath_slope = vector of bathymetric slopes corresponding to each shoreline point
%   shoreline = shoreline coordinates [x y]
% Outputs:
%   Xmesh = mesh in X-axis
%   Ymesh = mesh in Y-axis
%   Zmesh = mesh in depth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract shoreline coordinates
x = shoreline(:,1);
y = shoreline(:,2);

% Plot shoreline
figure;
plot(x, y);
hold on;

% Create a grid around the lake basin
start_x = min(x) - 1e3;
end_x = max(x) + 1e3;
x_frac_step = (end_x - start_x) / 100;

start_y = min(y) - 1e3;
end_y = max(y) + 1e3;
y_frac_step = (end_y - start_y) / 100;

x_extent = start_x:x_frac_step:end_x;
y_extent = start_y:y_frac_step:end_y;

[X, Y] = meshgrid(x_extent, y_extent);

% Identify points inside the lake boundary
polyin = polyshape([x y]);
polyout = polybuffer(polyin, 2e3);

[in, on] = inpolygon(X, Y, polyout.Vertices(:,1), polyout.Vertices(:,2));

plot(X(~in), Y(~in), 'ro');
plot(X(in | on), Y(in | on), 'go');

X_lake = X(in | on);
Y_lake = Y(in | on);

% Calculate bathymetry depths
min_dis = zeros(size(X_lake)); % Distance to the nearest shoreline point
assigned_slope = zeros(size(X_lake)); % Corresponding slope for each point

for j = 1:length(X_lake)
    % Compute distances to all shoreline points
    distances = sqrt((x - X_lake(j)).^2 + (y - Y_lake(j)).^2);
    [min_dis(j), imin] = min(distances);
    assigned_slope(j) = bath_slope(imin); % Use the slope of the closest shoreline point
end

% Calculate depths based on distances and assigned slopes
min_dep = assigned_slope .* min_dis;

% Interpolate depths to a regular grid
xLon = linspace(min(X_lake), max(X_lake), 1E+3);
yLat = linspace(min(Y_lake), max(Y_lake), 1E+3);
[Xmesh, Ymesh] = meshgrid(xLon, yLat);
zDep = griddata(X_lake, Y_lake, min_dep, Xmesh, Ymesh, 'linear');

% Assign NaN to points outside the lake basin
[in, on] = inpolygon(Xmesh, Ymesh, x, y);
zDep(~in & ~on) = NaN;

% Visualize results
figure;
mesh(Xmesh, Ymesh, zDep);
hold on;
[Mcon, Ccon] = contour3(Xmesh, Ymesh, zDep, [0 0], 'k', 'LineWidth', 2, 'ShowText', 1, 'LabelSpacing', 2000);
h = plot(x, y, '--k', 'LineWidth', 2);
z = get(h, 'ZData');
set(h, 'ZData', z + 10);
view(0, 90); % XY plane view
title('Bathtub Model for Lake Depth with Variable Slope');
colorbar;

end
