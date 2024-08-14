function [freqs,f_step] = make_frequency_vector(Model)
% MAKES A VECTOR OF LOG-SPACED FREQUENCIES BEING USED IN THE WAVE MODEL
% INPUT: 
%   Model  : containing information on min, max, and number of frequencies to include
% OUTPUT:
%   freqs  : frequency vector for spectrum
%   f_step : frequency step size of spectrum

    f_step=(log(Model.max_freq)-log(Model.min_freq))/(Model.Fdim-1);             
    freqs = exp(log(Model.min_freq)+(0:Model.max_freq-1)*f_step);  

end