function sum_Qs = calc_Qs_waves(my_shoreline,Wind)
% all input angles are radians

% returns Qs of waves for shoreline my_shoreline for wind climate Wind

rel_angle = deg2rad(-179:180);  

% Wind.phi0 = Wind direction
% Wind.E_pdf = fractional time of energy in phi0 direction

% Wind.H0 = wave height in wind direction phi0
Wind.H0 = 1;
% Wind.T0 = wave height in wind direction phi0
Wind.T0 = 1;

% Compute Wave Energy PDF and Sediment Transport Function
for m = 1:numel(my_shoreline)

    LST = CERC(Wind.H0, Wind.T0, rel_angle);  % LST from waves

    kappa = 10;
    % PDF rotated to POV of shoreline
    E_pdf = vonMises(Wind.phi0 - my_shoreline(m), kappa, rel_angle);  % Energy PDF (Von Mises)

    % Sediment transport (Qs)
    Qs = E_pdf.* LST;  % Combined sediment transport

    theta =  rel_angle + my_shoreline(m); % all possible angles given a starting shoreline orientation
    valid_theta = (theta >= -pi/2) & (theta <= pi/2);
    %Qs(~valid_theta) = NaN;

    % Find min and max Qs
    [minQ(m), iminQ(m)] = min(Qs);
    [maxQ(m), imaxQ(m)] = max(Qs);

    sum_Qs(m) = maxQ(m) - minQ(m);  % Net transport max

end



end

