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
time_to_run = 60*2;                                                         % time to run model
wind_direction = pi;                                                     % wind direction

[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.02);

zDep(zDep<1) = NaN;
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction(1),zDep,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);  

make_input_map(Planet,Model,Wind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRIDDED PLOT FOR DEPTH
x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);

num_cells_x = Model.LonDim;                                                % number of cells in x direction
num_cells_y = Model.LatDim;                                                % number of cells in y direction

WIDTH = (x_max-x_min)/num_cells_x;                                             % Width of each grid cell (# cells * width of 1 cell = width of all cells)
HEIGHT = (y_max-y_min)/num_cells_y;                                             % Height of each grid cell


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

% find theta
for i = 1:numel(x)
    if i < numel(x)
        shoreface_angle = atan2(y(i+1) - y(i), x(i+1) - x(i));
    else % wrap around to the first point to complete the circle
        shoreface_angle = atan2(y(1) - y(i), x(1) - x(i));
    end

    % Calculate the relative wave angle
    wave_front_angle = Wind.dir + pi / 2;
    theta(i) = mod(2*pi - (wave_front_angle - shoreface_angle), 2 * pi);

    % Limit theta to the range [0, pi/2]
    if theta(i) >= pi/2 || theta(i) <= 0
        theta(i) = NaN;
    end
end

% PLOT WAVE GRID ON SHOREFACE POINTS

figure;
plot(x,y,'-k')
hold on

for i = 1:num_cells_x
    for j = 1:num_cells_y

        % boundaries of grid cell
        x_lower = x_min + (i-1) * WIDTH;
        x_upper = x_min + i * WIDTH;
        y_lower = y_min + (j-1) * HEIGHT;
        y_upper = y_min + j * HEIGHT;

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
    cell_x = floor((x(i) - x_min) / WIDTH) + 1;
    cell_y = floor((y(i) - y_min) / HEIGHT) + 1;
    
    % skip points outside coordinate grid
    if cell_x < 1 || cell_x > num_cells_x || cell_y < 1 || cell_y > num_cells_y
        continue  
    end
    
    % get H, T of current cell
    H_cell = sig_wave.H(cell_y, cell_x);
    T_cell = sig_wave.T(cell_y, cell_x);
    
    if isnan(H_cell)
        [H_cell, T_cell] = check_neighbor(sig_wave, cell_x, cell_y, num_cells_x, num_cells_y, x(i), y(i), x_min, y_min, WIDTH, HEIGHT);
    end

    if isnan(H_cell)
        fprintf('no valid neighbors found %i\n',i)
    end
    % Calculate sediment transport rate Qs based on theta(i)
    Qs(i) = CERC(H_cell, T_cell, theta(i));
    psi(i) = shoreline_stability(H_cell, T_cell, theta(i),Model.bathy_map(cell_x,cell_y));
end

plot(x(~isnan(theta)),y(~isnan(theta)),'or','MarkerFaceColor','r')
plot(x(isnan(theta)),y(isnan(theta)),'ok','MarkerFaceColor','k')

scatter(x,y,100,Qs,"filled")
title('Qs')

% dQ_dx = NaN(size(Qs));
% for i = 1:numel(Qs)
% 
%     if isnan(Qs(i))
%         continue
%     end
% 
%     if i < numel(Qs)
%         dx = sqrt((y(i+1)-y(i))^2 + (x(i+1)-x(i))^2);
%         dQ_dx(i) = (Qs(i+1) - Qs(i))/dx;
% 
%     else
%         dx = sqrt((y(1)-y(i))^2 + (x(1)-x(i))^2);
%         dQ_dx(i) = (Qs(1) - Qs(i))/dx;
% 
%     end
% end
% 
% figure;
% scatter(x,y,100,dQ_dx,"filled")
% colorbar
% title('dQ/dx')


figure;
scatter(x,y,100,psi,"filled")
title('\Psi')
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

