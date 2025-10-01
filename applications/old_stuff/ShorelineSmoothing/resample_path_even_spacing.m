function [x_even, y_even] = resample_path_even_spacing(x, y, spacing)
%RESAMPLE_PATH_EVEN_SPACING Resample (x,y) path with even spacing
%
%   [x_even, y_even] = resample_path_even_spacing(x, y, spacing)
%
%   Inputs:
%     x, y     - Vectors of original coordinates (irregularly spaced)
%     spacing  - Desired distance between resampled points
%
%   Outputs:
%     x_even, y_even - Resampled coordinates with even spacing

    % Ensure column vectors
    x = x(:);
    y = y(:);

    % Compute cumulative arc length along the original path
    dx = diff(x);
    dy = diff(y);
    segment_lengths = hypot(dx, dy);
    arc_length = [0; cumsum(segment_lengths)];

    % Define new arc-length positions for even spacing
    total_length = arc_length(end);
    arc_length_even = (0:spacing:total_length)';

    % Interpolate x and y to the new arc-length positions
    x_even = interp1(arc_length, x, arc_length_even, 'linear');
    y_even = interp1(arc_length, y, arc_length_even, 'linear');
end
