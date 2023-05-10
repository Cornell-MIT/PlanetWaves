clc
clear
close all

tic
windspeeds = [13 20];

m = 31;                                                                    % number of gridpoints in x-direction
n = 15;                                                                    % number of gridpoints in y-direction

bathy_map = 100.*ones(m,n);

% from titanpool at 90K
% 
% rho_methane = 540;
% nu_methane = 3e-7; % m2/s
rho_water = 997;
nu_water = 1e-6;



[sigH_flat,T0] = makeWaves(windspeeds,rho_water,nu_water,bathy_map,1,100); % [m]

figure;
plot(1:length(sigH_flat),sigH_flat(1,:),'-o')
hold on;
plot(1:length(sigH_flat),sigH_flat(2,:),'-o')
hold off
xlabel('time step')
ylabel('sig H [m]')
legend('13 m/s','20 m/s')

for i = 1:length(T0)
    for j = 1:length(windspeeds)
        T0_time(j,i) = T0{j,i}(5,5);
    end
end

PM_T_13 = (2*pi)/((0.877*9.81)/13);
PM_T_20 = (2*pi)/((0.877*9.81)/20);


figure
plot(1:length(T0),T0_time(1,:),'-o')
hold on
plot(1:length(T0),T0_time(2,:),'-o')
plot(1:length(T0),PM_T_13.*ones(size(1:length(T0))),'--k')
plot(1:length(T0),PM_T_20.*ones(size(1:length(T0))),'--k')
xlabel('time step')
ylabel('T [sec]')
legend('13 m/s','20 m/s','PM 13 m/s','PM 20 m/s')
toc