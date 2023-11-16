clc
clear
close all

lm = readtable('TAM_Ligeia_timeseries_lat80.3_lon112.5.txt');

u = 0.76.*lm.u_m_s_;
v = 0.76.*lm.v_m_s_;

figure;
scatter(u,v)
hold on
plot(mean(u),mean(v),'o','MarkerFaceColor','r')
axis([-3 3 -3 3])
title('Ligeia Mare')

angle = rad2deg(atan2(v,u));
magU = sqrt(u.^2 + v.^2);

figure;
subplot(2,1,1)
histogram(angle)
xlabel('Direction [degree from East]')
title('Ligeia Mare')
subplot(2,1,2)
histogram(magU)
xlabel('|U|')