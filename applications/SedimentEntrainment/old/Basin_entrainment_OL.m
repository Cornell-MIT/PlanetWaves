clc
clear
close all

addpath(fullfile('..','..','data/Titan/TitanLakes/altimeter_pass'))
addpath(fullfile('..','..','data/Titan/TitanLakes/shoreline'))
addpath(fullfile('..','..','data/Titan/TitanLakes/Bathymetries/bathtub_bathy'))
A_slope = 2e-3; % section A from altimeter Hayes2010
L_slope = 1e-3; % section L from altimeter Hayes2010


sand_OL = 6.6160; % for 1.4 m/s winds

load('ontario_shoreline_pts.mat')
lon_shore = mod(ontarioshorelinepts.x,360);
lon_shore = [lon_shore ;lon_shore(1)];
lat_shore = ontarioshorelinepts.y;
lat_shore = [lat_shore ;lat_shore(1)];


[x_shore, y_shore, distances,ref_pt] = titan_geodesic_distances(lat_shore,lon_shore);

% [A_Xmesh,A_Ymesh,A_zDep] = make_bathtub_lake(A_slope,[x_shore y_shore]);
% [L_Xmesh,L_Ymesh,L_zDep] = make_bathtub_lake(L_slope,[x_shore y_shore]);
load('A_slope.mat','A_Xmesh','A_Ymesh','A_zDep')
load('L_slope.mat','L_Xmesh','L_Ymesh','L_zDep')

close all

load('T49_OL.mat')
x_path = mod(T49OL.x,360);
y_path = T49OL.y;

[x_49,y_49,~] = titan_geodesic_distances(y_path,x_path,ref_pt);

figure
contourf(A_Xmesh,A_Ymesh,A_zDep)
hold on
plot(x_49,y_49,'-r','LineWidth',2)
colorbar
title('A slope = 2e-3')

figure
contourf(L_Xmesh,L_Ymesh,L_zDep)
hold on
plot(x_49,y_49,'-r','LineWidth',2)
colorbar
title('L slope = 1e-3')

figure;
plot(x_shore,y_shore,'-k','LineWidth',2)
hold on
contour(A_Xmesh,A_Ymesh, A_zDep,[sand_OL sand_OL],'--k','LineWidth',2)
contour(L_Xmesh,L_Ymesh, L_zDep,[sand_OL sand_OL],'-k','LineWidth',2)
box on;
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
plot(x_49,y_49,'-r','LineWidth',2)


[A_x, A_y] = get_contour_xy(A_Xmesh,A_Ymesh,A_zDep,sand_OL);
[L_x, L_y] = get_contour_xy(L_Xmesh,L_Ymesh,L_zDep,sand_OL);

[A_lat,A_lon] = titan_xy_to_latlon(A_x{1},A_y{1},ref_pt);
[L_lat,L_lon] = titan_xy_to_latlon(L_x{1},L_y{1},ref_pt);

figure;
plot(lon_shore,lat_shore,'-k')
hold on
plot(A_lon,A_lat,'-k')
plot(L_lon,L_lat,'-k')

A_lon = fix_lon(A_lon);
A_basin = [A_lon' A_lat'];
L_lon = fix_lon(L_lon);
L_basin = [L_lon' L_lat'];

writematrix(A_basin,'A_basin.csv','Delimiter',',')
writematrix(L_basin,'L_basin.csv','Delimiter',',')


function new_lon = fix_lon(lon)

    for i = 1:numel(lon)
        if lon(i) > 180
            new_lon(i) = lon(i) - 360;
        else
            new_lon(i) = lon(i);
        end   
    end

end

function [x, y, distances, ref_point] = titan_geodesic_distances(lat, lon, ref_point)
%TITAN_GEODESIC_DISTANCES Computes UTM-like coordinates and distances on Titan
%   [x, y, distances, ref_point] = titan_geodesic_distances(lat, lon, ref_point)
%
%   Inputs:
%       lat, lon   - arrays of lat/lon in degrees
%       ref_point  - optional 2-element vector [lat0_rad, lon0_rad] (Titan radians)
%
%   Outputs:
%       x, y       - coordinates in meters relative to ref_point
%       distances  - cumulative great-circle distances in meters
%       ref_point  - [lat0_rad, lon0_rad], can be reused in future calls

    if numel(lat) ~= numel(lon)
        error('Latitude and longitude arrays must be the same length.');
    end
    if numel(lat) < 2
        x = 0; y = 0; distances = 0;
        if nargin < 3
            ref_point = [deg2rad(lat(1)), deg2rad(lon(1))];
        end
        return;
    end

    % Convert to radians
    lat_rad = deg2rad(lat(:));
    lon_rad = deg2rad(lon(:));
    
    % Titan's mean radius (m)
    R = 2574.36 * 1000;  % Zebker+2009

    % Handle optional reference point
    if nargin < 3 || isempty(ref_point)
        lat0 = lat_rad(1);
        lon0 = lon_rad(1);
        ref_point = [lat0, lon0];
    else
        lat0 = ref_point(1);
        lon0 = ref_point(2);
    end

    % Preallocate
    x = zeros(size(lat_rad));
    y = zeros(size(lat_rad));

    for i = 1:length(lat_rad)
        dlat = lat_rad(i) - lat0;
        dlon = lon_rad(i) - lon0;

        % Equirectangular projection
        x(i) = R * dlon * cos((lat_rad(i) + lat0)/2);
        y(i) = R * dlat;
    end

    % Haversine segment distances
    dlat = diff(lat_rad); dlon = diff(lon_rad);
    lat1 = lat_rad(1:end-1); lat2 = lat_rad(2:end);
    a = sin(dlat/2).^2 + cos(lat1).*cos(lat2).*sin(dlon/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    segment_distances = R * c;

    % Cumulative distance
    distances = [0; cumsum(segment_distances)];
end



function [lat_deg, lon_deg] = titan_xy_to_latlon(x, y, ref_point)
%TITAN_XY_TO_LATLON Converts local Titan XY back to lat/lon in degrees
%
%   [lat_deg, lon_deg] = titan_xy_to_latlon(x, y, ref_point)
%
%   Inputs:
%       x, y        - local coordinates in meters
%       ref_point   - [lat0_rad, lon0_rad] from titan_geodesic_distances
%
%   Outputs:
%       lat_deg, lon_deg - latitude and longitude in degrees (Titan)

    % Titan radius (m)
    R = 2574.36 * 1000;  % Zebker+2009

    lat0 = ref_point(1);
    lon0 = ref_point(2);

    % Convert back
    lat_rad = y / R + lat0;
    lon_rad = x ./ (R * cos((lat_rad + lat0)/2)) + lon0;

    % Convert to degrees
    lat_deg = rad2deg(lat_rad);
    lon_deg = rad2deg(lon_rad);
end

function [xc, yc] = get_contour_xy(X, Y, Z, level)
%GET_CONTOUR_XY Extracts (x, y) coordinates of contours at a given level
%   [xc, yc] = get_contour_xy(X, Y, Z, level)
%   X, Y: meshgrid coordinates
%   Z: scalar field (same size as X and Y)
%   level: the contour level you want to extract
%   Returns: xc, yc — cell arrays of x and y vectors for each contour line

    C = contourc(X(1,:), Y(:,1), Z, [level level]);

    xc = {};
    yc = {};
    k = 1;

    i = 1;
    while i < size(C, 2)
        this_level = C(1, i);
        n_points = C(2, i);

        if this_level == level
            x = C(1, i+1:i+n_points);
            y = C(2, i+1:i+n_points);
            xc{k} = x;
            yc{k} = y;
            k = k + 1;
        end

        i = i + n_points + 1;
    end
end

% x_OL = deg2utm_Titan(x_shore);
% % y_OL = deg2utm_Titan(y_shore);
% % 
% % figure;
% % plot(x_OL,y_OL)
% % shoreline(:,1) = x_OL;
% % shoreline(:,2) = y_OL;
% % [Xmesh,Ymesh,zDep] = make_bathtub_lake(rightmost_slope,shoreline);
% 
% 
% % figure;
% % plot(X_ol,Y_ol,'-k')
% % hold on
% % plot(X_cass,Y_cass,'or','MarkerFaceColor','r')
% % xlabel('Longitude ($^o$)','Interpreter','latex')
% % ylabel('Latitude ($^o$)','Interpreter','latex')
% 
% % ===== Grid Setup =====
% margin = 0.1;
% x_min = min([x_shore; x_path]) - margin;
% x_max = max([x_shore; x_path]) + margin;
% y_min = min([y_shore; y_path]) - margin;
% y_max = max([y_shore; y_path]) + margin;
% 
% gridRes = 1000;
% [xq, yq] = meshgrid(linspace(x_min, x_max, gridRes), linspace(y_min, y_max, gridRes));
% zq = nan(size(xq));
% 
% % ===== Lake Mask =====
% inLake = inpolygon(xq, yq, x_shore, y_shore);
% 
% % Determine which path points are inside the lake
% inLakePath = inpolygon(x_path, y_path, x_shore, y_shore);
% 
% entry_idx = find(diff(inLakePath) == 1) + 1; % enters lake
% exit_idx = find(diff(inLakePath) == -1);     % exits lake
% 
% shore_entry = [x_path(entry_idx), y_path(entry_idx)];
% shore_exit  = [x_path(exit_idx),  y_path(exit_idx)];
% 
% alt_res = titan_geodesic_distances(y_path, x_path);
% 
% 
% % Preallocate z_path
% z_path = nan(size(x_path));
% 
% for i = 1:length(x_path)
%     if i < entry_idx || i > exit_idx
%         continue; % skip points outside lake
%     else
%         z_path(i) = NaN;
%     end
% 
%     % distance from entry and exit along the path
%     d_entry = (i - entry_idx) * alt_res(i);
%     d_exit = (exit_idx - i) * alt_res(i);
% 
%     % sloping path in and out of basin
%     z_entry = leftmost_slope * d_entry;
%     z_exit = rightmost_slope * d_exit;
% 
%     % Choose minimum depth
%     z_path(i) = min(z_entry, z_exit);
% 
% 
% end
% 
% 
% for i = 1:numel(xq)
%     if ~inLake(i), continue; end
% 
%     % Current grid point
%     px = xq(i);
%     py = yq(i);
% 
%     % Compute squared distances to path points
%     d2 = (px - x_path).^2 + (py - y_path).^2;
% 
%     % Find nearest path point
%     [~, idx] = min(d2);
% 
%     % Assign depth
%     zq(i) = z_path(idx);
% end
% 
% % % ===== Plot =====
% figure;
% contourf(xq, yq, zq, 30, 'LineColor', 'none');
% hold on
% contour(xq,yq,zq,[10 10],'LineColor','r','LineWidth',2)
% plot(x_shore, y_shore, 'k-', 'LineWidth', 1.5);   % shoreline
% scatter(x_path, y_path, 10,z_path,'ko','filled')
% colorbar
% 
