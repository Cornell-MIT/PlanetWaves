function [m_f] = weighted_mean_freq(energy_spectrum,gridx,gridy,model)
% Computes the weighted mean frequency and corresponding direction
% INPUT: 
%   energy_spectrum = 4D energy spectrum with dimensions (Xgrid, Ygrid, Freq, Direction)
%   gridx           = location on grid in x dimension where frequency is evaluated
%   gridy           = location on grid in y dimension where frequency is evaluated
% OUTPUT:
%   m_f             = weighted mean frequency [Hz]

weighting_exp = 5; % 5 = M5, Reid1986; 8 = M8, Sobey+1986

directions = 0:(2*pi)/model.Dirdim:2*pi;
directions(end) = [];

f_vec = (log(model.max_freq) - log(model.min_freq)) / (model.Fdim - 1);             
f = exp(log(model.min_freq) + (0:model.Fdim-1) * f_vec);                        

% Extract the 2D frequency-direction spectrum at the specified grid point
freq_dir = squeeze(energy_spectrum(gridx, gridy, :, :));  % size: (Freq x Dir)

% Compute total energy for each frequency by integrating over direction
E_f = sum(freq_dir, 2);  % size: (Freq x 1)

% Compute the weighted mean frequency
m_f = sum(f(:) .* (E_f(:).^weighting_exp)) / sum(E_f(:).^weighting_exp);


end
