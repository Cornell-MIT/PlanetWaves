clc
clear
close all

EarthWavesTable = readtable('EarthData.csv');

xs = EarthWavesTable.WindSpeed_m_;
xs(isnan(xs)) = [];
[x,i] = sort(xs);
y = EarthWavesTable.WaveHeight_m_;
y(isnan(y)) = [];
y = y(i);
std_x = EarthWavesTable.STDWindSpeed_m_;
std_x(isnan(std_x)) = [];
sdt_x = std_x(i);
std_y = EarthWavesTable.STDWaveHeight_m_;
std_y(isnan(std_y)) = [];
std_y = std_y(i);

PM = (0.22.*x.^2)./9.81;

% Create scatter plot with error bars
figure;
scatter(x, y, 'filled');
hold on;
errorbar(x, y, std_y.*ones(size(y)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'r');
errorbar(x, y, std_x.*ones(size(x)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'b');
plot(x,PM,'-or','MarkerFaceColor','r')
grid on;

% Set plot title and labels
title('Earth Data');
xlabel('Wind Speed [m/s]');
ylabel('Wave Height [m]');
legend('Lake Superior (Deep)','Lake Superior (Deep) STD H','Lake Superior (Deep) STD |u|','Pierson-Moskowitz','location','best')