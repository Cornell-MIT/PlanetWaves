function grid_vals = smooth_lake_perimeter(grid_vals)
    % Get size of grid
    [rows, cols] = size(grid_vals);
    
    % Create a mask of valid (non-NaN) values
    valid_mask = ~isnan(grid_vals);
    
    % Find the boundary of the valid region
    boundary_mask = bwperim(valid_mask);  % Get perimeter of non-NaN region
    
    % Smooth only the boundary points
    for i = 2:rows-1
        for j = 2:cols-1
            if boundary_mask(i, j)  % If this point is on the perimeter
                % Get the valid neighbors
                neighbors = grid_vals(i-1:i+1, j-1:j+1);
                valid_neighbors = neighbors(~isnan(neighbors));
                
                % Replace with mean of valid neighbors
                if ~isempty(valid_neighbors)
                    grid_vals(i, j) = mean(valid_neighbors);
                end
            end
        end
    end
end
