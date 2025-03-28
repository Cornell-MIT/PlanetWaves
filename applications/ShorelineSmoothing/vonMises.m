function f = vonMises(mu,kappa,theta)
%   Inputs: 
%       mu: measure of location (analagous to mean in Gaussian), must be between 0:2pi
%       kappa: measure of concentration (1/k is analagous to variance in Gaussian), must be greater than zero
%           kappa = 1: uniform
%           kappa >> 1: highly concentrated aroud mu
%       theta: angle in polar coordinates
%   Outputs:
%       f: Von Mises probability density function for the angle mu

    f = exp(kappa*cos(theta-mu))/(2*pi*besseli(0,kappa));

end