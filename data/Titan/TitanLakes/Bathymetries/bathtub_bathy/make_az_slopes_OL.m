clc
clear
close all

load('OL_SHORELINE.mat','X_cor','Y_cor')

% TABLE 1 (MEAN DIP) Hahyes+2010 Bathymetry and Absorpitivity of Titan's Ontario Lacus
A_slope = 2e-3;
B_slope = 2.52e-3;
C_slope = 2.29e-3;
D_slope = 2.34e-3;
E_slope = 2.23e-3;
F_slope = 2.68e-3;
G_slope = 1.85e-3;
H_slope = 1.4e-3;
I_slope = 4.85e-3;
J_slope = 1.24e-3;
K_slope = 0.74e-3;
L_slope = 1e-3;
M_slope = 0.46e-3;

A = 2300;
B = 10;
C = 100;
D = 200;
E = 270;
F = 330;
G = 360;
H = 500;
I = 745;
J = 1180;
K = 1290;
L = 1575;
M = 1950;

x = X_cor;
x(1) = [];
y = Y_cor;
y(1) = [];
shoreline = [x y];

figure;
plot(x,y,'-k')
hold on
text(x(B),y(B),'B')
text(x(C),y(C),'C')
text(x(D),y(D),'D')
text(x(E),y(E),'E')
text(x(F),y(F),'F')
text(x(G),y(G),'G')
text(x(H),y(H),'H')
text(x(I),y(I),'I')
text(x(J),y(J),'J')
text(x(K),y(K),'K')
text(x(L),y(L),'L')
text(x(M),y(M),'M')
text(x(A),y(A),'A')


set(gca,'XDir','reverse')
set(gca,'YDir','reverse')

dips = [A A_slope; B B_slope; C C_slope; D D_slope; E E_slope; F F_slope; G G_slope; H H_slope; I I_slope; J J_soe]

% bath_slope1 = linspace(min_slope,max_slope,round(numel(x)/2));
% bath_slope2 = linspace(max_slope,min_slope,round(numel(x)/2));
% bath_slope = [bath_slope1 bath_slope2];

% [Xmesh, Ymesh, zDep] = make_bathtub_lake_w_azimuth_asymmetry(bath_slope, shoreline);
% 
% set(gca,'XDir','reverse')
% set(gca,'YDir','reverse')