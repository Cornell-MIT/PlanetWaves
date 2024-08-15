function [f,dlnf] = make_frequency_vector(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKES A VECTOR OF LOG-SPACED FREQUENCIES BEING USED IN THE WAVE MODEL
% INPUT: 
%   Model  : containing information on min, max, and number of frequencies to include
% OUTPUT:
%   freqs  : frequency vector for spectrum
%   f_step : frequency step size of spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dlnf=(log(model.max_freq)-log(model.min_freq))/(model.Fdim-1);           % frequency step size for log normal distribution
    f = exp(log(model.min_freq)+(0:model.Fdim-1)*dlnf);                      % frequencies for spectrum

end