clc
clear
close all

% make plots of signifigant wave height for slopes

windspeeds = [2 3];

m = 30;                                                                    % number of gridpoints in x-direction
n = 30;                                                                    % number of gridpoints in y-direction

ol_slope =10;

bathy_map = ones(m,n);
for i = 1:m
    bathy_map(i,:) = 100*(i);
end

surf(bathy_map)
% from titanpool at 90K

rho_liquid = 540;
nu_liquid = 3e-7; % m2/s
planet_gravity = 1.352;
planet_temp = 92;
planet_press = 1.5*101300;
surface_tension = 0.018;
gridX = 1000.0;
gridY = 1000.0;
time_step = 1;
num_time_steps = 100;
wind_dir = 0;


sigH_flat = makeWaves(windspeeds,wind_dir,rho_liquid,nu_liquid,planet_gravity,planet_temp,planet_press,surface_tension,bathy_map,gridX,gridY,time_step,num_time_steps); % [m]


figure;
plot(1:length(sigH_flat),sigH_flat(1,:),'-o')
hold on;
plot(1:length(sigH_flat),sigH_flat(2,:),'-o')
plot(1:length(sigH_flat),sigH_flat(3,:),'-o')
legend('u = 1 m/s','u = 2 m/s','u = 3 m/s')




