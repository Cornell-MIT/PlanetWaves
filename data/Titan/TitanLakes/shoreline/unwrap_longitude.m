function unwrapped_lon = unwrap_longitude(lon)
% avoid jumps at ±180° for longitude
%   Input:
%       lon - vector of longitude values (in degrees)
%   Output:
%       unwrapped_lon - adjusted longitudes for smooth plotting

    % makes column vector
    lon = lon(:);
    
    % Convert to radians for unwrapping
    lon_rad = deg2rad(lon);
    
    % Unwrap in radians
    unwrapped_rad = unwrap(lon_rad);
    
    % Convert back to degrees
    unwrapped_lon = rad2deg(unwrapped_rad);
end