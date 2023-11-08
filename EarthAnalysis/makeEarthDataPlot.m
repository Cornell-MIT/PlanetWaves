clc
clear
close all

EarthWavesTable = readtable('EarthData.csv');

x = EarthWavesTable.WindSpeed_m_;
y = EarthWavesTable.WaveHeight_m_;
std_x = EarthWavesTable.STDWindSpeed_m_;
std_y = EarthWavesTable.STDWaveHeight_m_;

PM = (0.22.*x.^2)./9.81;

% Create scatter plot with error bars
figure;
scatter(x, y, 'filled');
hold on;
errorbar(x, y, std_y.*ones(size(y)), 'vertical', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'r');
errorbar(x, y, std_x.*ones(size(x)), 'horizontal', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 10, 'Color', 'b');
plot(x,PM,'--r')
grid on;

% Set plot title and labels
title('Earth Data');
xlabel('Wind Speed [m/s]');
ylabel('Wave Height [m]');
legend('Lake Superior (Deep)','Pierson-Moskowitz','location','best')