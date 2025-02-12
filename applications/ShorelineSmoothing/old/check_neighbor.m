function [H_valid, T_valid, D_valid] = check_neighbor(sig_wave, cell_x, cell_y, num_cells_x, num_cells_y, x_point, y_point, x_min, y_min, width, height,Model,x,y)
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

    POI = 1;

    H_valid = NaN; T_valid = NaN; D_valid = NaN;
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
                    D_valid = Model.bathy_map(ny,nx);
                    min_dist = dist;  
                    if x_point == x(POI) && y_point == y(POI)
                        fprintf('x_center %i\n',x_center)
                        fprintf('y_center %i\n',y_center)
                        fprintf('D valid: %f\n',D_valid)
                    end
 
                end
            end
        end
    end
end



