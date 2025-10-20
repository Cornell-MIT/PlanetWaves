function fetch_results = calc_fetch(x, y, num_dir)
% calculate the fetch for each point along a closed shoreline
% Inputs:
% x, y: shoreline pixel coordinates for closed lake
% num_dir: number of directions (e.g., 360 for every degree)
%
% Output:
% fetch_results: matrix [num_points x num_dir] of fetch distances
% based on DDA algorithm originally implemented in python for Lake Superior

    % draw a buffer around the shoreline for the depth matrix
    buffer = 10;
    
    % make bounding box
    min_x = floor(min(x)) - buffer;
    max_x = ceil(max(x)) + buffer;
    min_y = floor(min(y)) - buffer;
    max_y = ceil(max(y)) + buffer;
    
    % Shift polygon coordinates so min is 1 for indexing
    x_shift = x - min_x + 1;
    y_shift = y - min_y + 1;
    
    % Calculate image size from bounding box
    rows = max_y - min_y + 1;
    cols = max_x - min_x + 1;
    
    % make binary depth mask from shifted polygon coordinates
    % liquid = 1, land = 0
    depth_mask = poly2mask(x_shift, y_shift, rows, cols);
    
    directions = linspace(0, 360 - 360/num_dir, num_dir);
    fetch_results = zeros(numel(x), num_dir);
    
    max_fetch = max(rows, cols) * 2; % maximum search length
    
    for i = 1:numel(x)
        start_x = x_shift(i);
        start_y = y_shift(i);
        
        for d = 1:num_dir
            deg = directions(d);
            dx = cosd(deg);
            dy = sind(deg);
            
            for step = 1:max_fetch
                x_test = round(start_x + dx * step);
                y_test = round(start_y + dy * step);
                
                % check if ray leaves the closed shoreline
                if x_test < 1 || x_test > cols || y_test < 1 || y_test > rows
                    fetch_results(i,d) = step;
                    break;
                end
                
                % Check if hit land (outside polygon)
                if depth_mask(y_test, x_test) == 0
                    fetch_results(i,d) = step;
                    break;
                end
            end
        end
    end
end
