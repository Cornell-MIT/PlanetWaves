clc
clear
close all

ATM_2_PASCAL = 101325;
RRR = 8.314; 

planet.rho_liquid = 588.15;                                            % TITANPOOL   
planet.nu_liquid = 8.084e-7;                                           % TITANPOOL   
planet.nua = 6.4e-9;                                                   % vapor N2 viscocity at 92 K (DIPPR)
planet.gravity = 1.352;                                                % Titan gravity
planet.surface_temp = 92;                                              % ~avg Titan surface temp
planet.surface_press = 1.5*ATM_2_PASCAL;                               % ~avg Titan surface pressure
planet.surface_tension = 0.032766;                                     % TITANPOOL
planet.kgmolwt = 0.028;                                                % (N2)

rho_a = planet.surface_press*planet.kgmolwt/(RRR*planet.surface_temp);
U10 = 0.4:0.1:4;
X = [0:10:180].*1000;
for i = 1:numel(U10)
    for j = 1:numel(X)
        H_sig(i,j) = (planet.nu_liquid^(0.002)/(planet.nua^(0.003)))*((rho_a/planet.rho_liquid)^0.36)*(((X(j)^0.07)*(U10(i)^1.6))/(planet.gravity^0.7));
    end
end

figure('Name','sig wave height from eqn')
contourf(X./1000,U10,H_sig);
colorbar
ylabel('wind speed m/s')
xlabel('Fetch km')

addpath(fullfile('..','..','\data\Titan\TitanLakes\Bathymetries\bathtub_bathy'))
addpath(fullfile('..','..','data/Titan/TAMwTopo/'))
addpath('C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\umwm_titan\applications\old_stuff\ShorelineSmoothing')
jet_wrap = vertcat(jet,flipud(jet)); % circular colormap

viridis_top = viridis;
viridis_top = viridis_top(1:end/2,:);
plasma_bottom = plasma;
plasma_bottom = plasma_bottom(end/2:end,:);
distinct_cmap = [viridis_top;plasma_bottom];

window_size_ang = 100;
Titan_cold.rho_s = 0.95*1000; 
Titan_cold.rho = 0.67*1000;
Titan_cold.nu = 0.003/10000;
Titan_cold.g = 1.35;

vidflumina.coldbedload_D50 = [3.8 10]/100;
vidflumina.susload_D50 = 6.35e-5;
vidflumina.width = [100 175];
vidflumina.slope = [0.0011 0.0015];

[river_susload,river_bedload] = riverine_flux(Titan_cold.rho_s,Titan_cold.rho,Titan_cold.nu,Titan_cold.g,vidflumina.coldbedload_D50(1),vidflumina.susload_D50,vidflumina.width(1),NaN);

load('ol_bathtub_0.46e_neg3.mat','X_cor','Y_cor','x_center','y_center','Xmesh','Ymesh','zDep')
X_cor(1) = [];
Y_cor(1) = [];
% shoreline is mapped CCW
x = (X_cor); 
y = (Y_cor);
[x,y] = make_circle(x,y);

figure('Name','window size for angle estimation');
plot(x,y,'-k','MarkerFaceColor','k')
hold on
scatter(x(1:window_size_ang),y(1:window_size_ang),50,jet(window_size_ang))
text(x(1)-1000,y(1),'start')
text(x(window_size_ang)-1000,y(window_size_ang),'end')
title('window size')


%load("shoreline_vs_Qs_wave.mat","my_shoreline","sum_Qs")
% OL_Qs_wave = interp1(wrap_180(rad2deg(my_shoreline)),sum_Qs,wrap_180(rad2deg(OL_shoreline_angle)));

Wind.phi0 = deg2rad(wrap_180(270)); % radians

OL_shoreline_angle = calc_shoreline_angle(x,y,window_size_ang); % radians

relative_angle = Wind.phi0 -  OL_shoreline_angle; % radians

OL_Qs_wave = calc_Qs_waves(OL_shoreline_angle,Wind); 

figure;
scatter(x,y,20,wrap_180(rad2deg(OL_shoreline_angle)),"filled")
hold on
quiver(0.6e5,-8e5,1e5*cos(deg2rad(Wind.phi0)),1e5*sin(deg2rad(Wind.phi0)))
colormap(jet_wrap)
title('shoreline angle')
colorbar


figure;
scatter(x,y,20,wrap_180(rad2deg(relative_angle)),"filled")
hold on
quiver(0.6e5,-8e5,1e5*cos(Wind.phi0),1e5*sin(Wind.phi0))
colormap(jet_wrap)
title('relative angle')
colorbar

figure;
scatter(x,y,20,OL_Qs_wave)
hold on
quiver(0.6e5,-8e5,1e5*cos(Wind.phi0),1e5*sin(Wind.phi0))
title('Qs wave')
colorbar








function [x_circle, y_circle] = make_circle(x,y)

    polyin = polyshape(x,y);
    [xc,yc] = centroid(polyin);

    for i = 1:numel(x)
        mydist(i) = sqrt((xc - x(i))^2 + (yc - y(i))^2);
    end

    my_r = max(mydist);

    theta = linspace(0, 2*pi, 1000);
    x_circle = xc + my_r.*cos(theta);
    y_circle = yc + my_r.*sin(theta);

end