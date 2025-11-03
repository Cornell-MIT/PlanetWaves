function Qs = calc_Diegaard_wave_flux(H0, L0, R, g, D50, nu, beta, alpha0)
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

    % Compute settling velocity using R
    w = calc_settling_velocity(R, g, D50, nu);
    % Dimensionless fall velocity
    w_star = w / sqrt(g * D50); % equation 13

    % Reference Shields parameter (Diegaarde 1986)
    theta_0 = 0.1 * (H0 / D50)^2.3 * sqrt(H0 / L0) * exp(-6.1 * w_star); % eqn 7

    Qs = NaN(size(alpha0));
    for a = 1:numel(alpha0)
        alpha_rad = alpha0(a);

        % Only compute if wave comes from seaward side
        if abs(alpha_rad) < pi/2

            alpha_norm = abs(alpha_rad) / (pi/2); % alpha0/90 deg in eqn 14

            % angular dependance for sed flux
            theta = ((sin(2.*alpha0.*(1-(0.4.*(alpha_norm).*(1-(alpha_norm)))))).^(5/2)).*theta_0; % eqn 14

            % Final sediment transport rate
            Qs(a) = theta * (H0 * sqrt(beta) * sqrt(R) * g * D50^3); % eqn 12
        end
    end

    function w = calc_settling_velocity(R, g, D50, nu)
        % CALC_SETTLING_VELOCITY
        %   Computes settling velocity using Dietrich (1982) 

        D_star = (R * g * D50^3) / (nu^2);
        log_D_star = log10(D_star);

        % Empirical fit for log10(W*)
        log_W_star = -3.76715 + ...
                      1.92944 * log_D_star - ...
                      0.09815 * log_D_star^2 - ...
                      0.00575 * log_D_star^3 + ...
                      0.00056 * log_D_star^4;

        W_star = 10^log_W_star;

        % Settling velocity
        w = (W_star * R * g * nu)^(1/3);
    end
end
