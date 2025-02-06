function dn_dt = shoreline_stability(H,T,theta,D)
% FIND SHORELINE STABILITY
    K2 = 1;
    dn_dt = ((K2/D)*(H^(12/5))*(T^(1/5)))*(cos(theta)^(1/5))*(cos(theta)^2 - (6/5)*sin(theta)^2);
    
end