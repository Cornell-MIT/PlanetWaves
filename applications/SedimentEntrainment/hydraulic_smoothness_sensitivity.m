clc
clear
close all


% VERTICAL VELOCITY PROFILE IN THE FLOW (affect of hydraulic roughness assumption)

% u(z) = (u_star/k)*ln(z/z0) <- Law of the Wall
%   z = from z0 to H
%   k = 0.4 (Von Karman Constant, dimensionless)
%   z0 = d_mean/30 <- assume rough turbulence (constant z0)

% CONSTANTS
d50 = [6.35e-5:1e-4:0.01];% diameters for [finegrain sand gravel] [m]
D_mean = mean(d50); %[m]
H = 15; % mean depth (good aprox for nearshore for all bath models in Vincent 2016)
k = 0.4;

% liquid
rho_w = 610; %kg/m3
kin_visc = 9e-7; %m2/s

% flows from prev models
u_coast =  0.046; %m/s 
u_strait = 0.64; %m/s 

% (1) At the nearshore, rough
z0_coast = D_mean/30; % rough turbulences
u_star_coast = (u_coast*k)/log(0.37*H/z0_coast); % using 4/10 rule to define u_star
z_coast = [z0_coast:16];
u_z_coast = (u_star_coast/k).*log(z_coast./z0_coast);
% (1b) At the nearshore, smooth
z0_coast_smooth = kin_visc/(9*u_star_coast); %smooth turbulence
u_star_coast_smooth = (u_coast*k)/log(0.37*H/z0_coast_smooth);
z_coast_smooth = [z0_coast_smooth:16];
u_z_coast_smooth = (u_star_coast_smooth/k).*log(z_coast_smooth./z0_coast_smooth);

% (2) At straits, rough
z0_strait = D_mean/30; %m, rough turbulences
u_star_strait = (u_strait*k)/log(0.37*H/z0_strait);
z_strait = [z0_strait:16];
u_z_strait = (u_star_strait/k).*log(z_strait./z0_strait);
% (2b) At straits, smooth
z0_strait_smooth = kin_visc/(9*u_star_strait); %smooth turbulence
u_star_strait_smooth = (u_strait*k)/log(0.37*H/z0_strait_smooth);
z_strait_smooth = [z0_strait_smooth:16];
u_z_strait_smooth = (u_star_strait_smooth/k).*log(z_strait_smooth./z0_strait_smooth);

fprintf('Average difference between rough and smooth for nearshore: %f\n',percent_difference(u_z_coast_smooth,u_z_coast))
fprintf('Average difference between rough and smooth for strait: %f\n',percent_difference(u_z_strait_smooth,u_z_strait))

% PLOTTING
figure('Name','z0 sensitivity')
subplot(1,2,1)
plot(u_z_coast,z_coast,'--r','LineWidth',3,'DisplayName','Rough')
hold on
plot(u_z_coast_smooth,z_coast_smooth,'-r','LineWidth',3,'DisplayName','Smooth')
xlabel('u (m/s)','fontsize',16)
ylabel('z (m)','fontsize',16)
grid on;
box on;
set(gca,'Xscale','log')
ax = gca;
ax.FontSize = 16;
legend('show','Location','best')
title('nearshore')
subplot(1,2,2)
plot(u_z_strait,z_strait,'--b','LineWidth',3,'DisplayName','Rough')
hold on;
plot(u_z_strait_smooth,z_strait_smooth,'-b','LineWidth',3,'DisplayName','Smooth')
xlabel('u (m/s)','fontsize',16)
ylabel('z (m)','fontsize',16)
grid on;
%legend('nearshore (0.46 m/s) rough','strait (0.64 m/s) rough','nearshore (0.46 m/s) smooth','strait (0.64 m/s) smooth')
set(gca,'Xscale','log')
ax = gca;
ax.FontSize = 16;
title('strait')
legend('show','Location','best')


function pd = percent_difference(x,y)


    pd = abs(x-y)./((x+y)./2);
    pd = pd.*100;
    pd = max(pd,[],'omitmissing');


end