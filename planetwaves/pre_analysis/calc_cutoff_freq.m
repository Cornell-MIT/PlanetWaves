function Model = calc_cutoff_freq(Planet,Model,Wind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATES A DYNAMIC CUTOFF FREQUENCY BETWEEN DIAGNOSTIC AND ADVECTING
% WAVELENGTHS THAT IS A FUNCTION OF THE PEAK PIERSON-MOSKOWITZ CURVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[freqs,~] = make_frequency_vector(Model);  
    for i = 1:length(Wind.speed)
        cutoff_value = (0.53*Planet.gravity)/Wind.speed(i);                    % calculate the cutoff freqency value [Hz]
        if Wind.speed(i) == 0                                                  % Check for division by zero
            Model.cutoff_freq(i) = 2;                                          % Assign 2 for undefined cutoff frequency
        else
            [~, index] = min(abs(freqs - cutoff_value));                       % Find the index of the closest frequency
            Model.cutoff_freq(i) = index;                                      % Store the index of the closest frequency
        end
    end

end