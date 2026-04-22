clc
clear
close all


addpath('C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\PlanetWaves\data\Titan\TitanLakes\shoreline')

load('ontariolacus_shoreline.mat','lon','lati')

x = lon(2:end);
y = lati(2:end);
Titan_Radius_m = 2575*1000;
[x, y] = deg2utm(y, x, Titan_Radius_m);

x = x(1:5:end);
y = y(1:5:end);

windows = [30];

for i = 1:numel(windows)
    windowsize = windows(i);
    sinuosity = calculate_sinuosity(x, y, windowsize);

    figure;
    scatter(x,y,50,sinuosity,"filled")
    title(["window size: ",num2str(windowsize)," points"])
    colorbar
end

failure = randn(length(x),1); % no correlation
negative_corr =  2 ./ sinuosity + 0.1*randn(length(x),1); % fake diffusivity that varies inverse to sin
positive_corr = sinuosity + 0.1*randn(length(x),1); % fake diffusivity that varies positive with sin

maxLag_test = min(50, floor(length(x)/10));

for test = 1:3
    if test == 1
        disp('This test should fail to find any correlation')
        diff = failure;
    elseif test == 2
        disp('This test should find a negative correlation')
        diff = negative_corr;
    elseif test == 3
        disp('This test should find a positive correlation')
        diff = positive_corr;
    end

    results = compare_diffusivity_sinuosity(x, y, diff, sinuosity, 'nperm',500,'maxLag', maxLag_test);
end