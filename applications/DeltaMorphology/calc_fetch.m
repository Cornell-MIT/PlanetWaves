function fetch_results = calc_fetch(x, y, POI,directions)
% Calculate the fetch for each point along a closed shoreline in specific directions.
% Inputs:
% x, y: shoreline pixel coordinates for closed lake
% POI : shoreline point finding fetches for along x and y
% directions: a scalar or vector of directions in degrees (e.g., 90 for east, or [0 90 180 270])
%
% Output:
% fetch_results: matrix [1 x num_directions] of fetch distances

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

    % Create binary depth mask from shifted polygon coordinates
    depth_mask = poly2mask(x_shift, y_shift, rows, cols);

    % Ensure directions is a row vector
    directions = directions(:)';
    num_dir = numel(directions);

    max_fetch = max(rows, cols) * 2; % maximum search distance

    for i = POI
        start_x = x_shift(i);
        start_y = y_shift(i);

        for d = 1:num_dir
            deg = directions(d);
            dx = cosd(deg);
            dy = sind(deg);

            for step = 1:max_fetch
                x_test = round(start_x + dx * step);
                y_test = round(start_y + dy * step);

                % If outside image bounds, stop
                if x_test < 1 || x_test > cols || y_test < 1 || y_test > rows
                    fetch_results(i, d) = step;
                    break;
                end

                % If hit land (outside polygon), stop
                if depth_mask(y_test, x_test) == 0
                    fetch_results(i, d) = step;
                    break;
                end
            end
        end
    end
end
