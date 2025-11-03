clc
clear
close all

make_plot = 0;

% Add path to shoreline data
addpath(fullfile('..','..','\data\Titan\TitanLakes\shoreline'))

% Load in shoreline
load('ontariolacus_shoreline.mat')
x = lon; y = lati;
x(1) = []; y(1) = [];

Titan_radius = 2574 * 1000; % meters
[x, y] = deg2utm(x, y, Titan_radius);

x = x ./ 1000;
y = y ./ 1000;

wind_direction = 270;

fetch_matrix = calc_fetch(x, y, 1:numel(x), wind_direction);

% index of points to plot with quiver
num_POI = 50;  % Number of points to sample
POI = randperm(numel(x), num_POI);  % Random unique indices

% Choose a colormap
cmap = jet(256);  % 256-color jet colormap

% Normalize fetch for color mapping
fetch_vals = fetch_matrix(POI, 1); % fetch along specified wind_direction
f_norm = (fetch_vals - min(fetch_vals)) / (max(fetch_vals) - min(fetch_vals)); % 0-1

figure; hold on; axis equal;
fill(x, y, [0.8 0.9 1], 'EdgeColor', 'b', 'LineWidth', 1.5); % Polygon
plot(x(POI), y(POI), 'ko', 'MarkerFaceColor', 'k');           % POIs

for i = 1:length(POI)
    start_x = x(POI(i));
    start_y = y(POI(i));
    len = fetch_matrix(POI(i), 1); % fetch along wind_direction
    angle = wind_direction;
    dx = len * cosd(angle);
    dy = len * sind(angle);

    % Map normalized fetch to colormap
    if len == 0
        col = [0.5 0.5 0.5];  % gray for zero fetch
    else
        idx = max(1, round(f_norm(i) * 255));  % map to 1-256
        col = cmap(idx, :);
    end

    % Draw arrow representing fetch
    quiver(start_x, start_y, dx, dy, 0, ...
        'Color', col, 'LineWidth', 2, 'MaxHeadSize', 0.5);

    % Annotate next to point
    text(start_x, start_y+1, sprintf('%.0f', len), ...
        'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end

xlabel('X ');
ylabel('Y');
grid on;
axis equal;
colorbar;  % optional, shows the colormap scale
