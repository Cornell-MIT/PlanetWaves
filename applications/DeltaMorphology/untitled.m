function [outputArg1,outputArg2] = Qs_waves(shoreline_angle,)
for m = 1:numel(my_shoreline)

    % Calculate relative wind direction
    wind = rel_angle - my_shoreline(m);  % Relative wind direction to shoreline

    % Filter valid wind directions within range [-pi/2, pi/2]
    valid_wind = (rad2deg(wind) > 90);

    % Von Mises PDF (you may need to define the vonMises function separately)
    E_pdf(m,:) = vonMises(phi0 - my_shoreline(m), kappa, wind);  % Energy PDF (Von Mises)

    % Compute CERC (sediment transport function, assuming it's a known function)
    LST = CERC(H0, T0, rel_angle);  % LST(wind-shoreline)

    % Sediment transport (Qs)
    Qs = E_pdf(m,:) .* LST;  % Combined sediment transport

    theta =  rel_angle + my_shoreline(m);
    valid_theta = (theta >= -pi/2) & (theta <= pi/2);
    Qs(~valid_theta) = NaN;

    % Find min and max Qs
    [minQ(m), iminQ(m)] = min(Qs);
    [maxQ(m), imaxQ(m)] = max(Qs);

    sum_Qs(m) = maxQ(m) - minQ(m);  % Net transport max

end
end

