function fetch = calc_fetch(x, y, wind_to_deg, min_fetch_distance)
% Calculates the fetch distance from all points on shoreline [x,y]
% along wind directions wind_to_deg
% Inputs:
%   x, y                : vectors of shoreline coordinates
%   wind_to_deg         : vector of wind directions in degrees (where wind blows to, degrees CCW from east)
%   min_fetch_distance  : minimum distance to consider a valid fetch (optional)
% Outputs:
%   fetch : [numel(wind_to_deg) x numel(x)] matrix of fetch distances
%
    make_plot = 0; 
    % Cyan: test point along ray
    % Red: invalid fetch
    % Yellow: ignored small segments
    % Green: valid fetch

    if nargin < 4
        min_fetch_distance = 1e-6; % minimum fetch distance that is considered non-zero 
    end

    

    num_points = numel(x); % number of shoreline points
    num_wind = numel(wind_to_deg); % number of wind directions
    fetch = NaN(num_wind, num_points);

    % Convert wind-to directions to wind-from directions (ray direction)
    wind_from_deg = mod(wind_to_deg + 180, 360); % wind coming from
    cos_wind = cosd(wind_from_deg);
    sin_wind = sind(wind_from_deg);

    % Threshold for small "leaks" out of the polygon that don't count
    max_leak_distance = 0.01 * range(x);
    ray_length = hypot(range(x), range(y)); % bounding box diagonal

    % make sure shoreline polygon is closed
    if x(1) ~= x(end) || y(1) ~= y(end)
        x(end+1) = x(1);
        y(end+1) = y(1);
    end


    % Loop over each shoreline point
    for pt_idx = 1:num_points
        start_point = [x(pt_idx); y(pt_idx)];

        if make_plot
            figure(1); clf;
            buffer = 0.05 * max([range(x), range(y)]);
            axis([min(x)-buffer, max(x)+buffer, min(y)-buffer, max(y)+buffer]);
            hold on;
            % fill in lake as grey
            fill(x, y, [0.9 0.9 0.9], 'EdgeColor', 'k', 'LineWidth', 1.2);
            xlabel('x'); ylabel('y');
            if pt_idx == 1
                gif('all_fetch.gif')
            end
            % Plot current shoreline point of interest as a blue dot
            plot(start_point(1), start_point(2), 'ob', 'MarkerFaceColor', 'b'); % start point
        end

        % Loop over each wind direction
        for wind_idx = 1:num_wind
            ray_vector = [cos_wind(wind_idx); sin_wind(wind_idx)];
            end_point = start_point + ray_length * ray_vector;

            % Small step along ray to check if initial direction is valid
            first_step = start_point + max_leak_distance * ray_vector;
            if make_plot
                % plot test point for determining if fetch immediatly leaves polygon
                plot(first_step(1), first_step(2), 'sc', 'MarkerFaceColor', 'c', 'MarkerSize', 8); % test point
            end

            [in, on] = inpolygon(first_step(1), first_step(2), x, y);
            if ~any(in) && ~any(on)
                % if test point is outside polygon then invalid direction
                fetch(wind_idx, pt_idx) = 0;
                if make_plot
                    % if test point is outside the polygon then fetch is invalid and cross out with a red x
                    plot(first_step(1), first_step(2), 'xr', 'MarkerSize', 10, 'LineWidth', 2); 
                    % plot invalid fetch as a red line
                    plot([start_point(1) end_point(1)], [start_point(2) end_point(2)], '-r', 'LineWidth', 1.2);
                end
                continue; % skip to next wind direction
            end

            % Find intersections with shoreline
            [xi, yi] = polyxpoly([start_point(1) end_point(1)], [start_point(2) end_point(2)], x, y);

            % Remove points too close to start (since the POI is also on the shoreline)
            if ~isempty(xi)
                dist_to_start = hypot(xi - start_point(1), yi - start_point(2));
                mask = dist_to_start >= 1e-10;
                xi = xi(mask); yi = yi(mask);
            end

            fetch(wind_idx, pt_idx) = 0; % default zero fetch

            if ~isempty(xi) % valid intersections found

                % Compute projected distance along ray
                distances = (xi - start_point(1))*ray_vector(1) + (yi - start_point(2))*ray_vector(2);
                [sorted_distances, sort_idx] = sort(distances);
                xi = xi(sort_idx); yi = yi(sort_idx);

                current_distance = 0;
                fetch_found = false;

                % Check segments between intersections
                for k = 1:length(sorted_distances)
                    segment_length = sorted_distances(k) - current_distance;

                    if segment_length < max_leak_distance
                        % if leak out of polygon is small, ignore
                        if make_plot
                            % plot distance being ignored as yellow line
                            plot([start_point(1)+current_distance*ray_vector(1), xi(k)],[start_point(2)+current_distance*ray_vector(2), yi(k)], '-y', 'LineWidth', 1.2);
                            % plot intersection for intersections being ignored as red circles
                            plot(xi(k), yi(k), 'or', 'MarkerFaceColor', 'r');
                        end
                        current_distance = sorted_distances(k);
                    else
                        % valid fetch
                        if sorted_distances(k) >= min_fetch_distance
                            fetch(wind_idx, pt_idx) = sorted_distances(k);
                            fetch_found = true;
                            if make_plot
                                % plot valid fetch as green line
                                plot([start_point(1), xi(k)], [start_point(2), yi(k)], '-g', 'LineWidth', 1.5);
                                % plot valid intersection as green circle
                                plot(xi(k), yi(k), 'og');
                            end
                        end
                        break;
                    end
                end

                if ~fetch_found && make_plot
                    % if no valid fetch was found, make all fetch red
                    plot([start_point(1) end_point(1)], [start_point(2) end_point(2)], '-r', 'LineWidth', 1.2);
                end
            else
                if make_plot
                    % make invalid fetches red
                    plot([start_point(1) end_point(1)], [start_point(2) end_point(2)], '-r', 'LineWidth', 1.2);
                end
            end
            if make_plot
                drawnow;
                gif;
            end

        end
        if make_plot
            pause(1) % pausing to make the gif easier to read
        end
    end
end
