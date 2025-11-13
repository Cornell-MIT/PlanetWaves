function [wave_height, wave_period, wave_length] = wind2wave(wind_speed, fetch, g)
% WIND2WAVE: Compute steady-state wave properties from wind & fetch
% and plot results cumulatively across multiple calls.
% wind_speed has size [1 x numel(timesteps)]
% fetch has size [numel(wind_speed) x numel(points)])
% returns waves with properties with size [numel(wind_speed) x numel(points)]

    persistent warning_shown 

    if isempty(warning_shown)
        warning('wind2wave: Using SMB/Young(1991) Fetch-Limited Growth with JONSWAP comparison.');
        warning_shown = true;
    end

    F = fetch;
    U = wind_speed(:);
    U = repmat(U, 1, size(F, 2)); % copies u where each row has same value of u to match fetch dimensions for vectorization

    % Dimensionless quantities
    X = (g .* F) ./ (U.^2); % fetch
    Hs_dimless = (5e-3) .* (tanh(0.0125 .* (X .^ 0.42))).^2; % wave height
    fp_dimless = 0.133 .* tanh(0.077 .* (X .^ 0.25)); % wave period

    % Physical quantities
    Hs = 4 .* sqrt(Hs_dimless) .* (U.^2 ./ g); 
    fp = fp_dimless .* (g ./ U);
    Tp = 1 ./ fp;
    Tp(isinf(Tp)) = NaN;
    Lp = (g .* Tp.^2) ./ (2 * pi);

    wave_height = Hs;
    wave_period = Tp;
    wave_length = Lp;

end
