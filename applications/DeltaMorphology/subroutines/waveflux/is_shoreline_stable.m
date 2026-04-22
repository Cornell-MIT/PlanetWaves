function dn_dt = is_shoreline_stable(H,T,ang)
% calculates the sum of shoreline stability (positive = bump is erased, negative = bump is exaggerated)  for a single spot on the shoreline
% INPUTS:
%   H = wave height time series (m) [1 x time]
%   T = wave period time series (s) [1 x time] 
%   ang = relative angle of approach (psi0 - theta) (deg) [1 x time]
% OUTPUTS:
%   dn_dt = summed normalized stability for a single spot on the shoreline with angle theta [1 x 1]

    make_plot = 0;

    % Convert to radians 
    ang_rad = deg2rad(ang);  

    % make all of them columns
    H = H(:); T = T(:); ang_rad = ang_rad(:);

   
    dn_dt= NaN(size(H));

    % shoreline instability (Ashton&Murray2006) -- equation 8
    valid_idx = ang_rad <= pi/2 & ang_rad >= 0;
    dn_dt(valid_idx) = angular_dependance(H(valid_idx),T(valid_idx),ang_rad(valid_idx));
    %reflecting function over zero
    reflect_idx = ang_rad >= -pi/2 & ang_rad < 0;
    dn_dt(reflect_idx) = angular_dependance(H(reflect_idx),T(reflect_idx),-ang_rad(reflect_idx));

    if make_plot
        figure('Name','Shoreline Stability');
        plot(ang,dn_dt)
        xlim([-180 180])
    end

    dn_dt = sum(dn_dt,'omitnan');


    function dn_dt = angular_dependance(H,T,ang_vals)
        % only computing relative value, so setting K2 and D to 1 for all points
        K2 = 1;
        D = 1;
        dn_dt = ((K2/D).*(H.^(12/5)).*(T.^(1/5))).*(cos(ang_vals).^(1/5)).*(cos(ang_vals).^2 - (6/5).*sin(ang_vals).^2);
    end

end