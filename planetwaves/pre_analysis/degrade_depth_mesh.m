function [depth_resized, buoy_loc_resized, new_resolution, Xmesh_resized, Ymesh_resized, newpath_x, newpath_y] = ...
    degrade_depth_mesh(depth, buoy_loc, grid_resolution, resizeFactor, Xmesh, Ymesh, path_x, path_y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to degrade resolution of a depth grid and related data
%   INPUT:
%       depth            : 2D array of depths [m]
%       buoy_loc         : grid location [longitude latitude] of buoy
%       grid_resolution  : [X-resolution Y-resolution] of depth grid [m]
%       resizeFactor     : multiplicative factor to degrade dimensions
%       Xmesh, Ymesh     : Original meshgrid coordinates
%       path_x, path_y   : Coordinates of a path
%   OUTPUT:
%       depth_resized    : 2D array of depths in degraded dimensions
%       buoy_loc_resized : Resized buoy location
%       new_resolution   : [X-resolution Y-resolution] of resized grid
%       Xmesh_resized    : X-coordinates for resized grid
%       Ymesh_resized    : Y-coordinates for resized grid
%       path_resized     : Resized path coordinates [x, y]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Handle NaN values in depth
    depth(isnan(depth)) = 0;
    
    % Calculate new grid resolution
    gridX_res_new = grid_resolution(1) * (1 / resizeFactor);
    gridY_res_new = grid_resolution(2) * (1 / resizeFactor);
    new_resolution = [gridX_res_new, gridY_res_new];
    
    % Resize depth array
    depth_resized = imresize(depth, resizeFactor, "bicubic");
    
    % Adjust buoy location
    buoy_loc_resized = buoy_loc .* resizeFactor;
    
    % Generate resized mesh grids
    original_x = Xmesh(1, :);
    resized_x = linspace(min(original_x), max(original_x), size(depth_resized, 2));
    original_y = Ymesh(:, 1);
    resized_y = linspace(min(original_y), max(original_y), size(depth_resized, 1));
    [Xmesh_resized, Ymesh_resized] = meshgrid(resized_x, resized_y);
    
    % Rescale path coordinates
    t = linspace(0, 1, numel(path_x)); % Parameterize original path
    new_t = linspace(0, 1, round(numel(path_x) * resizeFactor)); % Adjust for resolution
    newpath_x = interp1(t, path_x, new_t, 'linear');
    newpath_y = interp1(t, path_y, new_t, 'linear');

end
