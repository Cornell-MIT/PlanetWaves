clc
clear
close all

% making waves within a basin of uniform depth everywhere instead of sloping

% WAVES FOR 8 DIRECTIONS
wind_direction = 0:45:315;
wind_direction = deg2rad(wind_direction);
fne = {'0deg.mat','45deg.mat','90deg.mat','135deg.mat','180deg.mat','225deg.mat','270deg.mat','315deg.mat'};

load('asylake1.mat')
load('asylake1_bathtub.mat')

