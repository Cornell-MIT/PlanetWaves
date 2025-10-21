function [wave_height,wave_period,wave_length] = wind2wave(wind_speed,fetch)

% takes in properties of wind and computes steady state wave properties
% for volumetric flux calculations
warning('wind2wave: Using Fetch and Duration Limited Growth empirical functions, need to update to our model')

    g = 1.352;
    for i = 1:numel(wind_speed)
        u = wind_speed(i);
        for j = 1:numel(fetch)
            x = fetch(j);
            
            % dimensionless fetch
            X = (g*x)/(u^2);

            % SMB empricial functions
            calc_H_dimless = @(X) 0.30 .* (1 - exp(-3.3e-4 .* X.^0.5));
            calc_Tp_dimless = @(X) 5.0 .* (1 - exp(-1.0e-3 .* X.^0.3));
    
            Hs(i,j) = calc_H_dimless(X) .* u.^2 ./ g;
            Tp(i,j) = calc_Tp_dimless(X).* u./ g;
            Lp(i,j) = g .* Tp(i,j).^2 ./ (2 * pi);
    
        end
    end

    wave_height = Hs;
    wave_period = Tp;
    wave_length = Lp;
    
end