clc
clear
close all

% PLOT WAVES IN ONTARIO LACUS

addpath(fullfile('..','..','planetwaves'))  
addpath(fullfile('..','..','planetwaves','pre_analysis'))
addpath(fullfile('..','..','planetwaves','post_analysis'))

addpath(fullfile('..','..','data','Titan','TitanLakes','Bathymetries','bathtub_bathy'))
load('..\..\data\Titan\TitanLakes\Bathymetries\bathtub_bathy\ol_bathtub_0.002000_slope.mat','zDep')
load('..\..\data\Titan\TitanLakes\shoreline\OL_SHORELINE.mat','X_cor','Y_cor')
%load('..\data\Titan\TitanLakes\Bathymetries\bathtub_bathy\OL_shoreline_depth.mat','lon','lat','zDep')

zDep = imrotate(zDep,180);
% isolate main basin of interest
%zDep = smoothed_ol;
% zDep = imrotate(zDep,180);
% zDep = imrotate(zDep,-90);
% zDep(:,80:end) = [];
% zDep(1:95,:) = [];
% zDep(90:end,:) = [];
% zDep = imrotate(zDep,-90);
% 
% zDep_orig = zDep;


% MODEL INPUTS
planet_to_run = 'Titan-OntarioLacus';
buoy_loc = [630 255];                                                        % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
test_speeds = [3];                                                         % wind speed
time_to_run = 120;                                                         % time to run model
wind_direction = [pi];                                                   % wind direction

[zDep,buoy_loc,grid_resolution] = degrade_depth_resolution(zDep,buoy_loc,grid_resolution,0.02);

zDep(zDep<1) = NaN;
[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction(1),zDep,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);  

make_input_map(Planet,Model,Wind)


for i = 1:numel(test_speeds)

    for j = 1:numel(wind_direction)
        Wind.speed = test_speeds(i);
        Wind.dir = wind_direction(j);

        Model = calc_cutoff_freq(Planet,Model,Wind);

        [~, htgrid{i,j}, ~, ~ , ~ , ~, PeakWave{i,j}] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  


    end
end

wave_H_grid = htgrid{1,1}{1,end}';

sig_wave = PeakWave{1};

x = X_cor;
y = Y_cor;

x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);

num_cells_x = Model.LonDim;
num_cells_y = Model.LatDim; 


A = (x_max-x_min)/num_cells_x;                                             % Width of each grid cell (# cells * width of 1 cell = width of all cells)
B = (y_max-y_min)/num_cells_y;                                             % Height of each grid cell


figure;
plot(x, y, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 1);  

shoreface_angle = NaN(size(x));
for i = 1:numel(x)-1
    shoreface_angle(i) = atan2(y(i+1) - y(i), x(i+1) - x(i));
end


wave_front_angle = Wind.dir + (pi)/2;
theta = wave_front_angle - shoreface_angle;
theta = mod(2*pi-theta, 2*pi);

% limit theta between 0 and pi/2
for i = 1:num_cells_x
    if theta(i) > pi/2 || theta(i) < 0
        theta(i) = NaN;
    end
end



num_colors = 100;
color1 = [229, 211, 82] / 255;
color2 = [172, 57, 49] / 255;

mycmap = [linspace(color1(1), color2(1), num_colors)', ...
                       linspace(color1(2), color2(2), num_colors)', ...
                       linspace(color1(3), color2(3), num_colors)'];


figure;
scatter(x,y,10,rad2deg(theta),"filled")
hold on

point_indices = NaN(size(x));
for i = 1:num_cells_x
    for j = 1:num_cells_y

        % boundaries of grid cell
        x_lower = x_min + (i-1) * A;
        x_upper = x_min + i * A;
        y_lower = y_min + (j-1) * B;
        y_upper = y_min + j * B;
        
        %maxH = max(max(wave_H_grid));
        %H_color = round((wave_H_grid(j,i)/maxH)*100);
        maxH = max(max(Model.bathy_map));
        H_color = round((Model.bathy_map(j,i)/maxH)*100);
  
        if isnan(H_color)
            H_color = 1;
            rr = rectangle('Position', [x_lower, y_lower, A, B], 'EdgeColor', 'k','FaceColor',mycmap(H_color,:),'FaceAlpha',0); 

        end
        % draw grid
        rr = rectangle('Position', [x_lower, y_lower, A, B], 'EdgeColor', 'k','FaceColor',mycmap(H_color,:),'FaceAlpha',0.5); 
        
        % find all points within grid
        in_cell = (x >= x_lower & x < x_upper) & (y >= y_lower & y < y_upper);
        
        % index based on the grid cell (i, j)
        index = (i-1) * num_cells_y + j;
        point_indices(in_cell) = index;
        
        [qs_row,qs_col] = ind2sub(size(Model.bathy_map),index);

        if any(in_cell)

            % all points witin a single cell will have the same deepwater wave properties but can vary in angle
            theta_cell = shoreface_angle(in_cell);

            H_cell = sig_wave.H(qs_row,qs_col);
            T_cell = sig_wave.T(qs_row,qs_col);

            for pp = 1:numel(theta_cell)
                % sediment transport (based on CERC)s
                Qs(pp) = ((H_cell).^(12/5)).*((T_cell).^(1/5)).*(cos(theta_cell(pp)).^(6/5)).*sin(theta_cell(pp));
            end
        end
    end
end


xlim([x_min - A, x_max + A]);
ylim([y_min - B, y_max + B]);
hold off




dQ_dx = zeros(size(Qs));

for i = 1:numel(Qs)-1
    new_dif = Qs(i+1) - Qs(i);
    dQ_dx(i) = sum([dQ_dx(i) new_dif],"omitmissing");
    dQ_dx(i+1) = sum([dQ_dx(i) new_dif],"omitmissing");
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_sorted, y_sorted] = orderbytheta(x, y)
    % orders points (x, y) by their angle theta relative to center

    center_x = mean(x);
    center_y = mean(y);

    angles = atan2(y - center_y, x - center_x);
    [~, sort_idx] = sort(angles); 
   
    x_sorted = x(sort_idx);
    y_sorted = y(sort_idx);

end

