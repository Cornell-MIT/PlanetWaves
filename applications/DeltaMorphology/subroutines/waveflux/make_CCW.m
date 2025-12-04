function [x_CCW,y_CCW] = untitled(x,y)
% If the shoreline points are clockwise (CW), reverses the point order

    % Compute polygon signed area to detect orientation
    signed_area = 0.5 * sum(x .* circshift(y, -1) - y .* circshift(x, -1));
    
    if signed_area < 0
            % Polygon is clockwise -> reverse order to make it CCW
            x_CCW = flip(x);
            y_CCW = flip(y);
        else
            % Already CCW
            x_CCW = x;
            y_CCW = y;
    end

end

