function [Eastward, Northward] = deg2utm(lat_deg, lon_deg, planet_radius, lon0_deg)
% DEG2UTM Convert lat/lon (deg) to UTM-like coordinates for spherical planet
% Inputs:
%   lat_deg       - Latitude(s) in degrees
%   lon_deg       - Longitude(s) in degrees (can be any range)
%   planet_radius - Radius of the planet (meters)
%   lon0_deg      - (Optional) Central meridian in degrees (default = auto zone-based)
%
% Outputs:
%   Eastward  - UTM-like Easting (meters)
%   Northward - UTM-like Northing (meters)

    % Normalize longitude to [-180, 180)
    lon_deg = mod(lon_deg + 180, 360) - 180;

    % central meridian as the mean of unwrapped longitudes
    lon_rad_all = deg2rad(lon_deg);
    lon0_rad = angle(mean(exp(1i * lon_rad_all)));
    lon0_deg = rad2deg(lon0_rad);

    % Convert degrees to radians
    lat_rad = deg2rad(lat_deg);
    lon_rad = deg2rad(lon_deg);
    lon0_rad = deg2rad(lon0_deg);

    % UTM scale factor guess this will work for now
    k0 = 0.9996;

    % spherical Transverse Mercator projection
    B = cos(lat_rad) .* sin(lon_rad - lon0_rad);
    x = 0.5 * planet_radius * log((1 + B) ./ (1 - B));
    y = planet_radius * atan(tan(lat_rad) ./ cos(lon_rad - lon0_rad));

    Eastward = k0 * x + 500000;
    Northward = k0 * y;

    % False northing for southern hemisphere
    Northward(lat_deg < 0) = Northward(lat_deg < 0) + 10000000;
end
