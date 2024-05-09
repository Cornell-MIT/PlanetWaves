function [depth_resized,buoy_loc_resized,new_resolution] = degrade_depth_resolution(depth,buoy_loc,grid_resolution,resizeFactor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function will take a depth profile and resize it to a lower resolution 
%   INPUT:
%       depth            : 2D array of depths [m]
%       buoy_loc         : grid location [longitude latitude] of buoy
%       grid_resolution  : [X-resolution Y-resolution] of depth grid (aka pixel width and height) [m]
%       resizeFactor     : multiplicative factor to degrade dimensions by
%  OUTPUT:
%       depth_resized    : 2D array of depths in degraded dimensions
%       buoy_loc_resized : grid location of buoy in reshaped array
%       new_resolution   : [X-resolution Y-resolution] in reshaped array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    depth(isnan(depth)) = 0;
    gridX_res_new = (grid_resolution(1))*(1/resizeFactor);
    gridY_res_new = (grid_resolution(2))*(1/resizeFactor);
    new_resolution = [gridX_res_new gridY_res_new];
    
    pos = [buoy_loc(1), buoy_loc(2)];
    depth = imresize(depth, resizeFactor, "bicubic");
    depth_resized = round(depth);
    
    buoy_loc_resized = ceil(pos * resizeFactor);


end