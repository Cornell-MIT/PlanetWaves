clc
clear
close all

% testing functionality of makeWaves function

windspeeds = 0.4:1:3.3;

m = 31;                                                                    % number of gridpoints in x-direction
n = 15;                                                                    % number of gridpoints in y-direction
    
rho_liquid = 465;
nu_liquid = 0.0031/1e4;

bathy_map = 100.*ones(m,n);

sigH = makeWaves(windspeeds,rho_liquid,nu_liquid,bathy_map);