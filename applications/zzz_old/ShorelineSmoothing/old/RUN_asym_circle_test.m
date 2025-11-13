clc
clear
close all
% 
% % SHORELINE SMOOTHING ON ASYMETRIC CIRCLE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % HOUSEKEEPING
addpath(fullfile('..','..','planetwaves'))  
addpath(fullfile('..','..','planetwaves','pre_analysis'))
addpath(fullfile('..','..','planetwaves','post_analysis'))


% simple bathtub model of depths
[x,y,Xmesh, Ymesh, Depth] = make_fake_shoreline(2.4852e+05,2000,2e4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP BATHYMETRY FOR WAVE MODEL


% MODEL INPUTS
planet_to_run = 'Titan-OntarioLacus';
buoy_loc = [4500 500];                                                      % grid location [x,y]
grid_resolution = [1000 1000];                                             % pixel width and pixel height [m]
test_speeds = 3;                                                           % wind speed
time_to_run = 60;                                                          % time to run model
wind_direction = 0;                                                        % wind direction

figure;
subplot(1, 2, 1);
surf(Xmesh, Ymesh, Depth,'EdgeColor','none');
view(2)
hold on;
plot3(x, y, max(Depth(:)) * ones(size(x)), 'r', 'LineWidth', 2);
title('Original Depth and Path');
xlabel('X'); ylabel('Y'); zlabel('Depth');

[zDep, buoy_loc, new_resolution, Xmesh, Ymesh, ~, ~] = degrade_depth_mesh(Depth, buoy_loc, grid_resolution, 1/50, Xmesh, Ymesh, x, y);

% Display results
% zDep(zDep<0) = NaN;

buoy_loc = [10 10];
zDep = round(zDep);
zDep(zDep==0) = NaN;


[Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction(1),zDep,buoy_loc);

Model.gridX = grid_resolution(1);                                              
Model.gridY = grid_resolution(2);  

subplot(1, 2, 2);
pcolor(Xmesh, Ymesh, zDep);
colorbar
hold on;
plot3(x, y, max(zDep(:)) * ones(size(x)), 'r', 'LineWidth', 2);
title('Resized Depth and Path');
xlabel('X'); ylabel('Y'); zlabel('Depth');
view(2)
colorbar


figure
pcolor(Xmesh,Ymesh,zDep)
view(2)
title('original')
colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRIDDED PLOT FOR DEPTH

WIDTH = (max(x)-min(x))/Model.LonDim;                                      % Width of each grid cell (# cells * width of 1 cell = width of all cells)
HEIGHT = (max(y)-min(y))/Model.LatDim;                                     % Height of each grid cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN WAVE MODEL (need H, T) OF MAIN ENERGY COMPONENT
Wind.speed = test_speeds;
Wind.dir = wind_direction;

Model = calc_cutoff_freq(Planet,Model,Wind);

make_input_map(Planet,Model,Wind)

[~,~,~,~,~,~,PeakWave] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  

% EXTRACT GRID OF WAVE HEIGHTS AND PERIODS
sig_wave = invert_attribute(PeakWave);


POI = 1;
h = plot_wave_grid(x,y,WIDTH,HEIGHT,Model,sig_wave);
hold on
plot(x(POI),y(POI),'.r')

sig_wave.H = smooth_lake_perimeter(sig_wave.H);
sig_wave.T = smooth_lake_perimeter(sig_wave.T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE Qs and Psi for each sub-grid point from wave grid
[wind_dir,wave_front_angle,shoreline_angle,relative_angle] = calc_shoreline_angle(x,y,Wind);

wave_energy = calc_wave_energy(x,y,sig_wave,WIDTH,HEIGHT,Model,Planet,relative_angle);
[Qs,mydiff] = shoreline_smoothing(x,y,sig_wave,relative_angle,WIDTH,HEIGHT,Model);
sinu = calc_sinuosity(x,y,200);


figure;
scatter(x,y,100,sinu,"filled")
colorbar
title('sinuosity')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE POINTS NORMAL RELATIVE TO THE DIRECTION OF THE WIND

my_angle = point_rel2_wind(x,y,wind_dir);

% jet_wrap = vertcat(jet,flipud(jet));
% figure;
% 
% scatter(x,y,100,rad2deg(my_angle),"filled")
% colormap(jet_wrap);
% colorbar




alongwind_diff = mydiff(along_wind);
alongwind_sin = sinu(along_wind);
not_alongwind_diff = mydiff(~along_wind);
not_alongwind_sin = sinu(~along_wind);
alongwind_wave_energy = wave_energy(along_wind);
not_alongwind_wave_energy = wave_energy(~along_wind);

r_along_diffsin = r_squared(alongwind_diff,alongwind_sin,'alongwind diff','alongwind sin');
r_notalong_diffsin = r_squared(not_alongwind_diff,not_alongwind_sin,'not alongwind diff','not alongwind sin');
r_along_energysin = r_squared(alongwind_wave_energy,alongwind_sin,'alongwind wave energy','alongwind sin');
r_notalong_energysin = r_squared(not_alongwind_wave_energy,not_alongwind_sin,'not alongwind wave energy','alongwind sin');

figure;
scatter(x,y,100,mydiff,'ok')
hold on
scatter(x(along_wind),y(along_wind),100,alongwind_diff,"filled")
colorbar
title('diffusivity')

[~,dans_along,xans_along,yans_along] = kde2d([alongwind_diff' alongwind_sin']);

[~,dans_notalong,xans_notalong,yans_notalong] = kde2d([not_alongwind_diff' not_alongwind_sin']);
figure;
contour3(xans_notalong,yans_notalong,dans_notalong,50,'-r')
view(2)
hold on;
contour3(xans_along,yans_along,dans_along,50,'-k')

xlabel('diffusivity')
ylabel('sinuoisty')
ylim([0 1])



figure;
scatter(x,y,100,wave_energy,'ok')
hold on
scatter(x(along_wind),y(along_wind),100,alongwind_wave_energy./max(alongwind_wave_energy),"filled")
colorbar
title('wave energy')

[~,dans_along,xans_along,yans_along] = kde2d([alongwind_wave_energy' alongwind_sin']);

[~,dans_notalong,xans_notalong,yans_notalong] = kde2d([not_alongwind_wave_energy' not_alongwind_sin']);
figure;
contour3(xans_notalong,yans_notalong,dans_notalong,50,'-r')
view(2)
hold on;
contour3(xans_along,yans_along,dans_along,50,'-k')

xlabel('wave energy')
ylabel('sinuoisty')
ylim([0 1])

function STRUCT = invert_attribute(STRUCT)

fields = fieldnames(STRUCT);
for i = 1:numel(fields)
    if isnumeric(STRUCT.(fields{i})) && ismatrix(STRUCT.(fields{i}))
        STRUCT.(fields{i}) = transpose(STRUCT.(fields{i}));
    end
end

end

function newval = closest_power2(val)
    newval = pow2(floor(log2(val)));
end

function [xEven, yEven] = even_spacing(x, y,gap)

    
    total_length = arclength(x,y,'linear');

    p = interparc(0:(gap/total_length):1,x,y,'linear');
    xEven = p(:,1);
    yEven = p(:,2);
    
end


function my_angle = point_rel2_wind(x,y,wind_dir)
    polyin = polyshape(x,y);
    [xcent,ycent] = centroid(polyin);
    v1 = [cos(wind_dir), sin(wind_dir)];
    
    for pt = 1:numel(x)
        v2 = [(x(pt) - xcent), (y(pt) - ycent)];
        my_angle(pt) = angle_between(v1,v2);
    end

    figure;
    plot(x,y)
    hold on
    quiver(xcent,ycent,2e4*v1(1),2e4*v1(2),'-r');
    plot([xcent x(pt)],[ycent y(pt)],'--g')

end


function myangle = angle_between(vec1,vec2)
    
    % normalize the vectors
    v1 = vec1/norm(vec1);
    v2 = vec2/norm(vec2);

    % find dot product
    dot_product = dot(v1,v2);
    % find cross product
    cross_product = det([v1; v2]);

    % find angle between difference vector
    myangle = atan2(cross_product, dot_product);

    % wrap between 0 and 2pi
    if myangle < 0
        myangle = myangle + 2*pi;
    end

end

