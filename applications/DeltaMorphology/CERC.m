function LST = CERC(H,T,ang)
% takes in angles in degree

ang_rad = deg2rad(ang);
valid_idx = (ang_rad >= -pi/2) & (ang_rad <= pi/2);

LST = NaN(size(ang));
LST(valid_idx) = (H^(12/5)) * (T^(1/5)) * ...
                 (cosd(ang(valid_idx)).^(6/5)) .* sind(ang(valid_idx));

end

