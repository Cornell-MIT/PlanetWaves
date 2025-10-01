clc
clear
close all

POI = 210;%50; %10

load('Model.mat')


% figure
% surf(Model.bathy_map)
% view(2)
% colorbar

% load('wind_0.mat')
% clear sin
% 
% figure;
% surf(imrotate(sig_wave.H,-90))
% set(gca,'XDir','reverse')
% view(2)
% title('H')
% colorbar
% 
% figure;
% surf(imrotate(sig_wave.T,-90))
% set(gca,'XDir','reverse')
% view(2)
% title('T')
% colorbar


load('..\..\data\Titan\TitanLakes\shoreline\OL_SHORELINE.mat','X_cor','Y_cor')
x = X_cor;y = Y_cor;
[x,y] = smooth_path(x,y,10);


% figure;
% scatter(x(~isnan(theta)),y(~isnan(theta)),100,psi(~isnan(theta)),'filled')
% hold on
% plot(x(POI),y(POI),'or')
% title('diffuse')
% axis([-5e4 20e4 -8.8e5 -6.8e5])
% colorbar;
% 
% figure;
% scatter(x,y,100,rad2deg(theta),'filled')
% hold on
% plot(x(POI),y(POI),'or')
% title('theta')
% axis([-5e4 20e4 -8.8e5 -6.8e5])
% colorbar

WIDTH = (max(x)-min(x))/Model.LonDim;                                     
HEIGHT = (max(y)-min(y))/Model.LatDim;  


plot_wave_grid(x,y,WIDTH,HEIGHT,Model)
hold on
plot(x(POI),y(POI),'or','MarkerFaceColor','r')

% NORTH
Wind.dir = 3*pi/2;
[wind_dir,wave_front_angle,shoreline_angle,relative_angle] = calc_shoreline_angle(x,y,Wind);

fprintf('wind: %0.1f deg\n',rad2deg(mod(pi/2 - Wind.dir, 2*pi)))
U = 1e4*cos(wind_dir);
V = 1e4*sin(wind_dir);
quiver(x(POI),y(POI),U,V,'k')

fprintf('wave: %0.1f deg\n',rad2deg(wave_front_angle(POI)))
U = 1e4*cos(wave_front_angle(POI));
V = 1e4*sin(wave_front_angle(POI));
quiver(x(POI),y(POI),U,V,'b')

fprintf('theta: %0.1f deg\n',rad2deg(shoreline_angle(POI)))
U = 1e4*cos(shoreline_angle(POI));
V = 1e4*sin(shoreline_angle(POI));
quiver(x(POI),y(POI),U,V,'g')

fprintf('relative angle: %0.1f deg\n',rad2deg(relative_angle(POI)))
U = 1e4*cos(relative_angle(POI));
V = 1e4*sin(relative_angle(POI));
quiver(x(POI),y(POI),U,V,'r')


% load("wind_pi.mat","sig_wave")
% 
% [Qs,psi] = shoreline_smoothing(x,y,sig_wave,relative_angle,WIDTH,HEIGHT,Model);
% 
% 
% figure;
% scatter(x,y,50,Qs,"filled")
% colorbar
% title('Qs')
% cmap = redblue(numel(psi));
% 
% figure;
% scatter(x,y,50,'k')
% hold on;
% scatter(x,y,50,-psi,'filled')
% 
% 
% clim([-2 2])
% colormap(cmap);
% colorbar
% title('diffusivity')
