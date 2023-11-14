clc
clear
close all

% [S,A] = shaperead('us_medium_shoreline.shp');
% 
% for i = 1:length(S)
%     a = S(i,:);
%     plot(a.X,a.Y)
%     title(i)
%     pause(1) % up to 1670
% end

lake = imread('superior_lld.tif');
lake(lake==-9999) = NaN;
lake(lake>=0) = NaN;