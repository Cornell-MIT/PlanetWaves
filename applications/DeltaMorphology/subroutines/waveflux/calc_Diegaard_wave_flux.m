function Qs = calc_Diegaard_wave_flux(H0, L0, rho_s,rho, g, D50, nu, beta, alpha0)
% CALC_DIEGAARD_WAVE_FLUX
%   Computes wave-driven sediment transport rate (Qs) based on the 
%   Diegaarde (1986) model.
%
% INPUTS:
%   H0      - Deep-water wave height (m)
%   L0      - Deep-water wave length (m)
%   R       - Relative submerged density (R = rho_s/rho - 1)
%   g       - Gravitational acceleration (m/s^2)
%   D50     - Median grain size (m)
%   nu      - Kinematic viscosity of water (m^2/s)
%   beta   - Beach slope (dimensionless)
%   alpha0  - Wave angle of approach relative to shore normal (radians)
%
% OUTPUT:
%   Qs      - Sediment transport rate for each angle in alpha0 (same size as alpha0)

    H0(isnan(H0)) = 0;
    L0(isnan(L0)) = 0;

    alpha0 = rad2deg(alpha0);

    % Compute settling velocity
    w = calc_settling_velocity(rho_s,rho, g, D50, nu);
    
    % Dimensionless fall velocity according to Diegaard
    w_star = w / sqrt(g * D50); % equation 13

    % net wave properties (Neinhuis 2015 supplement)
    Hnet = (sum(H0.^(2.8))/length(H0)).^(1/2.8);
    L0net = (sum(L0.^(0.5))/length(L0)).^(1/0.5);
    
    % Reference Shields parameter (Diegaarde 1986)
    theta_0 = 0.1 * ((Hnet / D50)^2.3) * sqrt(Hnet / L0net) * exp(-6.1 * w_star); % at alpha = 45; eqn 7

    Qs = NaN(size(alpha0));
    for a = 1:numel(alpha0)

        % Only compute if wave comes from seaward side
        alpha = alpha0(a);
    
        if abs(alpha) < 90 && alpha ~= 0
            
            alpha_norm = abs(alpha) / 90; 
        
            % angular dependence for sed flux
            sin_val = sind(2*abs(alpha)*(1 - ((0.4*alpha_norm)*(1 - alpha_norm))));
            
            theta_mag = (abs(sin_val).^(5/2)) * theta_0; % magnitude only
            
            theta = sign(alpha) * theta_mag; % reflect over zero to be CERC-like
        
            s = rho_s / rho; % relative density
            Qs(a) = theta * (Hnet * sqrt(beta) * sqrt((s-1) * g * (D50^3))); % eqn 12
        
        else
            Qs(a) = 0;
        end
        
    end

    function w = calc_settling_velocity(rho_s, rho, g, D50, nu)
        % CALC_SETTLING_VELOCITY
        %   Computes settling velocity using Dietrich (1982) 

        D_star = ((rho_s - rho)*g*(D50^3))/(rho*(nu^2));
        log_D_star = log(D_star);

        % Empirical fit for log10(W*)
        log_W_star = -3.76715 + ...
                      1.92944 * log_D_star - ...
                      0.09815 * log_D_star^2 - ...
                      0.00575 * log_D_star^3 + ...
                      0.00056 * log_D_star^4;

        W_star = exp(log_W_star);

        % Settling velocity
        w = (W_star * ((rho_s/rho)-1) * g * nu)^(1/3);
    end
end
