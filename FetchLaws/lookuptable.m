clc
clear
close all

tic
windspeeds = [13];
wind_dir = 0;
m = 31;                                                                    % number of gridpoints in x-direction
n = 15;                                                                    % number of gridpoints in y-direction

bathy_map = 100.*ones(m,n);

% from titanpool at 90K
% 
% rho_methane = 540;
% nu_methane = 3e-7; % m2/s
rho_water = 997;
nu_water = 1e-6;



sigH_flat = makeWaves(windspeeds,wind_dir,rho_water,nu_water,bathy_map,1,100); % [m]






toc