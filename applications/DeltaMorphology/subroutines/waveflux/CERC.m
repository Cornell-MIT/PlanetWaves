function Qs = CERC(H,T,rho,rho_s,g,ang)
% Calculates longshore sediment transport (LST) based on CERC formula

% ang: angles in degrees
    
    H(isnan(H)) = 0;
    T(isnan(T)) = 0;
    Qs = zeros(size(ang));

    % Convert to radians 
    ang_rad = deg2rad(ang);  
    valid_idx = (ang_rad >= -pi/2) & (ang_rad <= pi/2);
    
    % empirical constant
    K = 0.46*rho*(g^(3/2)); % Komar 1998
    p = 0.1; % dry mass void fraction (what to set this to?)
    gamma_b = 0.78; % ratio of breaking wave height to water depth
    n = 0.5; % for deepwater Cg = 0.5*C, 1 for shallow water
    K1 = (5.3e-6)*(K)*((1/(2*n))^(6/5))*(((sqrt(g*gamma_b))/(2*pi))^(1/5));

    % net wave properties (Neinhuis 2015 supplement)
    Hnet = (sum(H.^(12/5))/length(H)).^(5/12);
    Tnet = (sum(T.^(1/5))./length(T)).^(5);

    % CERC
    Qs(valid_idx) = K1*rho_s*(1-p)*(Hnet^(12/5))*(Tnet^(1/5)).*(cos(ang_rad(valid_idx)).^(6/5)) .* sin(ang_rad(valid_idx)); % coefficient * angular component

end
