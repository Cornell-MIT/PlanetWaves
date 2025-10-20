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

num_directions = 360;
fetch_matrix = calc_fetch(x, y, num_directions);
max_fetch = max(fetch_matrix,[],2,'omitmissing');

my_dir = linspace(0, 360 - 360/num_directions, num_directions);

figure;
imagesc(1:num_directions,my_dir,fetch_matrix)
xlabel('shoreline point')
ylabel('fetch km')

figure
plot(1:numel(max_fetch),max_fetch)
yline(12)
yline(62)
xlabel('shoreline point')
ylabel('max fetch km')

u_v_crit_fetch = [2 36*1000;
                  3 50*1000;
                  4 62*1000];

if make_plot
    figure('Name', 'Fetch');
    axis equal;
    hold on;
    
    for point_idx = 1:100:numel(x)
        x0 = x(point_idx);
        y0 = y(point_idx);
    
        % Clear axes for new plot
        cla;
    
        plot(x, y, 'k-', 'LineWidth', 2);  
        plot(x0, y0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');  
    
        angles = linspace(0, 2*pi, num_directions + 1);
        angles(end) = [];
    
        for j = 1:num_directions
            theta = angles(j);
            dist = fetch_matrix(point_idx, j);
            if ~isnan(dist)
                x_end = x0 + dist * cos(theta);
                y_end = y0 + dist * sin(theta);
                plot([x0, x_end], [y0, y_end], 'r--');
                text(x_end, y_end, num2str(dist),'Color','r');
            end
        end
    
        title(sprintf('Point %d of %d', point_idx, numel(x)));
        legend({'Shoreline', 'Start Point', 'Fetch', 'Intersection'}, 'Location', 'northeast');
        drawnow;
    
        % if point_idx == 1
        %     gif('ontario_lacus_fetch.gif','DelayTime',1);
        % else
        %     gif;
        % end
    end
end