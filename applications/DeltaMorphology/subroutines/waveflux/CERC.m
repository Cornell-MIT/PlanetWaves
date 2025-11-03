function LST = CERC(H,T,ang)
% Calculates longshore sediment transport (LST) based on CERC formula
% ang: angles in degrees

% Convert to radians 
ang_rad = deg2rad(ang);  
valid_idx = (ang_rad >= -pi/2) & (ang_rad <= pi/2);

LST = NaN(size(ang));

% CERC
LST(valid_idx) = (H(valid_idx).^(12/5)) .* (T(valid_idx).^(1/5)) ...
                 .* (cos(ang_rad(valid_idx)).^(6/5)) .* sin(ang_rad(valid_idx));
end
