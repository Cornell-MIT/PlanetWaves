clc
clear
close all

make_plot = 0;

% Add path to shoreline data
addpath(fullfile('..','..','\data\Titan\TitanLakes\shoreline'))

% Load in shoreline
load('ontariolacus_shoreline.mat')
x = lon; y = lati;
x(1) = []; y(1) = [];

% Calculate the regional shoreline orientaiton at each point
window_size_ang = 5; % size of window to define the regional shoreline orientation over
shoreline_angle = calc_regional_shoreline_angle(x, y,window_size_ang);
shoreline_angle = shoreline_angle';

% Juan's GCM winds
fl = 'C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\PlanetWaves\data\Titan\TAMwTopo';
load('OL_winds.mat','mag_wind','angle_wind')


Titan_radius = 2574 * 1000; % meters
[x, y] = deg2utm(x, y, Titan_radius);

x = x ./ 1000;
y = y ./ 1000;

num_directions = 360;

[E_pdf,theta] = make_wind_rose(mag_wind,angle_wind,1);

fetch_matrix = calc_fetch(x, y, num_directions);

crit_fetch = 36*1000;
max_height = calc_height(crit_fetch);


for point_idx = 1:numel(x)

    x0 = x(point_idx);
    y0 = y(point_idx);
    all_fetch = fetch_matrix(point_idx,:);
    wave_height = calc_height(all_fetch);
    wave_height(all_fetch>crit_fetch) = max_height;
    
    directional_energy = E_pdf .* wave_height;
    total_energy(point_idx) = sum(directional_energy);

    wave_sed_flux{point_idx} = calc_wave_flux(directional_energy,rad2deg(theta),rad2deg(shoreline_angle(point_idx)));

end

wave_sed_flux = cell2mat(wave_sed_flux);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate max flux from largest unobservable rivers

Titan_cold.rho_s = 0.95*1000; 
Titan_cold.rho = 0.67*1000;
Titan_cold.nu = 0.003/10000;
Titan_cold.g = 1.35;

vidflumina(1).width = 100;
vidflumina(1).slope = 0.0011;
vidflumina(2).width = 175;
vidflumina(2).slope = 0.0015;

for i = 1:2
    [river_susload(i),river_bedload(i)] = calc_riverine_flux(Titan_cold,vidflumina(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate fluvial dominance R (Qriver/Qwaves)

R = river_bedload(1).Qs./wave_sed_flux;

viridis_top = viridis;
viridis_top = viridis_top(1:end/2,:);
plasma_bottom = plasma;
plasma_bottom = plasma_bottom(end/2:end,:);
distinct_cmap = [viridis_top;plasma_bottom];

% Plot showing fluvial dominance ratio R 
% Large R -> River Dominated
% Small R -> Wave Dominated
% R = 1 -> Waves and rivers equal
figure('Name','Delta Morph R')
scatter(x,y,10,log10(R))
colormap(distinct_cmap)
colorbar;
clim([-2 2])
title('log10 R')
axis equal padded

function height = calc_height(fetch)

    height = fetch.^0.11;

end


    % tiledlayout(1,2)
    % 
    % nexttile
    % polarplot(deg2rad(theta), directional_energy,'-k');
    % hold on
    % polarplot(deg2rad(theta), E_pdf,'--r');
    % polarplot(deg2rad(theta),all_fetch./(max(all_fetch)*100),'--b')
    % title(sprintf('Directional Wave Energy at POI %s',num2str(point_idx)));
    % legend('Epdf * wave height','Epdf','fetch')
    % 
    % nexttile
    % plot(x,y,'-k')
    % hold on
    % plot(x(point_idx),y(point_idx),'or','MarkerFaceColor','r')
    % plot(x, y, 'k-', 'LineWidth', 2);  
    % plot(x0, y0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');  
    % 
    % angles = linspace(0, 2*pi, num_directions + 1);
    % angles(end) = [];
    % 
    % for j = 1:num_directions
    %     theta = angles(j);
    %     dist = fetch_matrix(point_idx, j);
    %     if ~isnan(dist)
    %         x_end = x0 + dist * cos(theta);
    %         y_end = y0 + dist * sin(theta);
    %         plot([x0, x_end], [y0, y_end], 'r--');
    %         text(x_end, y_end, num2str(dist),'Color','r');
    %     end
    % end