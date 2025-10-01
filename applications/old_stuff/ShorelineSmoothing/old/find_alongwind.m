function along_wind = find_alongwind(x,angle_from_shore)
for i = 1:numel(x)
    if angle_from_shore(i) > pi/2 && angle_from_shore(i) < 3*pi/2
        along_wind(i) = false;
    else
        along_wind(i) = true;
    end
end
end