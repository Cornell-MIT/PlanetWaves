clc
clear
close all

% SHORELINE SMOOTHING AT ONTARIO LACUS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOUSEKEEPING
addpath(fullfile('..','..','planetwaves'))  
addpath(fullfile('..','..','planetwaves','pre_analysis'))
addpath(fullfile('..','..','planetwaves','post_analysis'))
addpath(fullfile('..','..','data','Titan','TitanLakes','Bathymetries','bathtub_bathy'))
% simple bathtub model of depths
load('..\..\data\Titan\TitanLakes\Bathymetries\bathtub_bathy\ol_bathtub_0.002000_slope.mat','zDep')
% hi-res coordinates of shoreline
load('..\..\data\Titan\TitanLakes\shoreline\OL_SHORELINE.mat','X_cor','Y_cor')
x = X_cor;y = Y_cor;
[x,y] = smooth_path(x,y,10);
%[x,y] = make_circle(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURES STUFF
num_colors = 100;
color1 = [229, 211, 82] / 255;
color2 = [172, 57, 49] / 255;

mycmap = [linspace(color1(1), color2(1), num_colors)', ...
                       linspace(color1(2), color2(2), num_colors)', ...
                       linspace(color1(3), color2(3), num_colors)'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP BATHYMETRY FOR WAVE MODEL

zDep = imrotate(zDep,180);

% MODEL INPUTS
planet_to_run = 'Titan-OntarioLacus';
buoy_loc = [630 255];                                                      % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
test_speeds = 3;                                                         % wind speed
time_to_run = 60*10;                                                         % time to run model
wind_direction = pi;                                                     % wind direction

[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.06);

zDep(zDep<1) = NaN;
[a,b] = size(zDep);
% zDep = max(max(zDep)).*ones(20,20);

[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction(1),zDep,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);  

make_input_map(Planet,Model,Wind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRIDDED PLOT FOR DEPTH

WIDTH = (max(x)-min(x))/Model.LonDim;                                             % Width of each grid cell (# cells * width of 1 cell = width of all cells)
HEIGHT = (max(y)-min(y))/Model.LatDim;                                             % Height of each grid cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN WAVE MODEL (need H, T) OF MAIN ENERGY COMPONENT

for i = 1:numel(test_speeds)

    for j = 1:numel(wind_direction)
        Wind.speed = test_speeds(i);
        Wind.dir = wind_direction(j);

        Model = calc_cutoff_freq(Planet,Model,Wind);

        [~, ~, ~, ~ , ~ , ~, PeakWave{i,j}] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  


    end
end

% EXTRACT GRID OF WAVE HEIGHTS AND PERIODS
test_u_ind = 1;
test_dir_ind = 1;
sig_wave = PeakWave{test_u_ind,test_dir_ind}';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE Qs for each sub-grid point from wave grid

% Initialize arrays for theta and Qs, matching the order of shoreline points (x, y)
theta = NaN(size(x));  % Relative angle of attack
Qs = NaN(size(x));     % Sediment transport rate
psi = NaN(size(x));
% find theta
for i = 1:numel(x)
    if i < numel(x)
        shoreface_angle(i) = atan2(y(i+1) - y(i), x(i+1) - x(i));
    else % wrap around to the first point to complete the circle
        shoreface_angle(i) = atan2(y(1) - y(i), x(1) - x(i));
    end
    
    shoreface_angle(i) = wrapToPi(shoreface_angle(i));
    % Calculate the relative wave angle
    wave_front_angle = wrapToPi(Wind.dir + pi / 2);
    theta(i) = wrapToPi(wave_front_angle - shoreface_angle(i));

    % % Limit theta to the range [0, pi/2]
    if abs(theta(i)) > pi/2 || abs(theta(i)) < 0
        theta(i) = NaN;
    end
end

% PLOT WAVE GRID ON SHOREFACE POINTS

figure;
plot(x,y,'-k')
hold on

for i = 1:Model.LonDim
    for j = 1:Model.LatDim

        % boundaries of grid cell
        x_lower = min(x) + (i-1) * WIDTH;
        x_upper = min(x) + i * WIDTH;
        y_lower = min(y) + (j-1) * HEIGHT;
        y_upper = min(y) + j * HEIGHT;

        % PLOT WAVE GRID ON SUB-GRID ATTACK ANGLE
        % using relative size compared to max for coloring grid
        maxH = max(max(Model.bathy_map));
        H_color = round((Model.bathy_map(j,i)/maxH)*100);
        
        % draw grid
        if isnan(H_color)
            rr = rectangle('Position', [x_lower, y_lower, WIDTH, HEIGHT],'FaceColor',[1 1 1],'EdgeColor', 'k','FaceAlpha',0); 
        else
            rr = rectangle('Position', [x_lower, y_lower, WIDTH, HEIGHT], 'EdgeColor', 'k','FaceColor',mycmap(H_color,:),'FaceAlpha',0.5); 
        end
    end
end

% Calculate Qs and diffusivity for each (x, y) point based on the same theta ordering
for i = 1:numel(x)

    % skip if POI is NaN
    if isnan(x(i))
        continue
    end
    
    % find grid containing POI (x(i), y(i)) 
    cell_x = floor((x(i) - min(x)) / WIDTH) + 1;
    cell_y = floor((y(i) - min(y)) / HEIGHT) + 1;
    
    % skip points outside coordinate grid
    if cell_x < 1 || cell_x > Model.LonDim || cell_y < 1 || cell_y > Model.LatDim
        %continue  
        if cell_y > Model.LatDim
            cell_y = Model.LatDim;
        elseif cell_x > Model.LonDim
            cell_x = Model.LonDim;
        else
            continue
        end
    end
    
    % get H, T of current cell
    H_cell = sig_wave.H(cell_y, cell_x);
    T_cell = sig_wave.T(cell_y, cell_x);
    
    if isnan(H_cell)
        [H_cell, T_cell] = check_neighbor(sig_wave, cell_x, cell_y, Model.LonDim, Model.LatDim, x(i), y(i), min(x), min(y), WIDTH, HEIGHT);
    end

    if isnan(H_cell)
        fprintf('no valid neighbors found %i\n',i)
    end
    % Calculate sediment transport rate Qs based on theta(i)
    Qs(i) = CERC(H_cell, T_cell, theta(i));
    psi(i) = shoreline_stability(H_cell, T_cell, theta(i),zDep(cell_x,cell_y));
end

plot(x(~isnan(theta)),y(~isnan(theta)),'or','MarkerFaceColor','r')
plot(x(isnan(theta)),y(isnan(theta)),'ok','MarkerFaceColor','k')

scatter(x,y,100,Qs,"filled")
colorbar
title('Qs')


figure;
scatter(x,y,100,psi,"filled")
colorbar
title('\Psi')

sin = calc_sinuosity(x,y,5);

figure;
scatter(x,y,100,sin,"filled")
colorbar
title('sinuosity')

figure;
plot(psi,sin,'ok','MarkerFaceColor','k')
xlabel('\Psi')
ylabel('sin')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Qs = CERC(H,T,theta)
   % FIND SEDIMENT FLUX FOR WAVES OF HEIGHT (H), PERIOD (T), WITH ATTACK ANGLE (theta)
   % INPUT:
   %    H = wave height [m]
   %    T = wave period [s]
   %    theta = wavefront angle - shoreface angle 
   % OUTPUT:
   %    Qs = Volumetric Sediment Flux

    Qs = ((H).^(12/5)).*((T).^(1/5)).*(cos(theta).^(6/5)).*sin(theta);
end

function dn_dt = shoreline_stability(H,T,theta,D)
% FIND SHORELINE STABILITY
    K2 = 1;
    dn_dt = -((K2/D)*(H^(12/5))*(T^(1/5)));
    dn_dt = dn_dt*(cos(theta)^(1/5))*(cos(theta)^2 - (6/5)*sin(theta)^2);
    
end
function [H_valid, T_valid] = check_neighbor(sig_wave, cell_x, cell_y, num_cells_x, num_cells_y, x_point, y_point, x_min, y_min, width, height)
    % CHECK EIGHT NEAREST NEIGHBORS IF NO WAVES PRESENT EXACTLY IN GRID CELL WILL RETURN NEAREST non-NAN 
    % INPUT:
    %   sig_wave    : grid of wave heights
    %   cell_x      : current cell row
    %   cell_y      : current cell column
    %   num_cells_x : number of cells in x direction
    %   num_cells_y : number of cells in y direction
    %   x_point     : POI x-value
    %   y_point     : POI y-value
    %   x_min       : minimum values of x in overall shoreline
    %   y_min       : minimum values of y in overall shoreline
    %   width       : width of grid cell
    %   height      : height of grid cell
    % OUTPUT:
    %   H_valid     : nearest-neighbor non-NaN waveheight
    %   T_valid     : nearest-neighbor non-NaN wave period

    H_valid = NaN; T_valid = NaN;
    min_dist = inf;  
    
    % Eight nearest neighbors
    neighbor_offsets = [-1, -1; -1, 0; -1, 1;
                        0, -1;        0, 1;
                        1, -1;  1, 0; 1, 1];
    
    
    for k = 1:size(neighbor_offsets, 1)
        nx = cell_x + neighbor_offsets(k, 1);
        ny = cell_y + neighbor_offsets(k, 2);
        
        % check bounds
        if nx >= 1 && nx <= num_cells_x && ny >= 1 && ny <= num_cells_y
            % Check if non-NaN H 
            if ~isnan(sig_wave.H(ny, nx))
                
                x_center = x_min + (nx - 0.5) * width;
                y_center = y_min + (ny - 0.5) * height;
                
                % find distance POI to neighbor's center
                dist = sqrt((x_point - x_center)^2 + (y_point - y_center)^2);
                
                % check if this cell is closer than prev valid cell and update if it is
                if dist < min_dist
                    H_valid = sig_wave.H(ny, nx);
                    T_valid = sig_wave.T(ny, nx);
                    min_dist = dist;  
                end
            end
        end
    end
end

function [xnew,ynew] = smooth_path(x,y,n)

numPoints = floor(length(x) / n);

xnew = zeros(1, numPoints);
ynew = zeros(1, numPoints);

for i = 1:numPoints
    startIdx = (i - 1) * n + 1;
    endIdx = startIdx + n - 1;
    xnew(i) = mean(x(startIdx:endIdx),"omitmissing");
    ynew(i) = mean(y(startIdx:endIdx),"omitmissing");
end

if xnew(1) ~= xnew(end) || ynew(1) ~= ynew(end)
    xnew = [xnew xnew(1)];
    ynew = [ynew ynew(1)];
end
end

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

end

function sinuosity = calc_sinuosity(x,y,ws)

    % distance walking along the shore
    for i = 1:length(x)
        
        if i + 1 <= length(x)
            each_diff(i) = sqrt((y(i+1) - y(i))^2 + (x(i+1)-x(i))^2);
        else 
            each_diff(i) = sqrt((y(1) - y(i))^2 + (x(1)-x(i))^2);
        end
    end


   for i = 1:length(x)
    
        if i - ws >= 1 && i + ws <= length(x)
            j1 = i - ws;
            j2 = i + ws;
            full_length_alongshore(i) = sum(each_diff(j1:j2));
        elseif i - ws >= 1 && i + ws > length(x)
            j1 = i - ws;
            j2 = i + ws - length(x);
            full_length_alongshore(i) = sum(each_diff(j1:length(x))) + sum(each_diff(1:j2));      
        elseif i - ws < 1 && i + ws <= length(x)
            j1 = length(x) + (i - ws);
            j2 = i + ws;
            full_length_alongshore(i) = sum(each_diff(j1:length(x))) + sum(each_diff(1:j2));
        else
            disp('missed points')
            disp('i')
        end
        
        % P1 - > P2 is distance as crow flies 
        pt1 = [x(j1),y(j1)]; 
        pt2 = [x(j2),y(j2)];
        
        abs_diff(i) = sqrt((pt2(2) - pt1(2))^2 + (pt2(1) - pt1(1))^2); % absolute distance between two points on the shoreline          
           
   end

   % sinuosity is between 0 and 1
    %   s -> 1 means less sinous
    %   s -> 0 means more sinous
    sinuosity = abs_diff./full_length_alongshore;
    

end