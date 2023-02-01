clc
clear 
close all

load('T_500.mat','UUvec','ht') % ran for Time = 500

figure
plot(UUvec(1:end-1),ht(1:end-2,end-1).*100,'-o')
hold on
xlabel('wind speed [m/s]')
ylabel('sig H [cm]')
grid on

load('T_10000.mat','UUvec','ht') % ran for Time = 10000


plot(UUvec(1:end-1),ht(1:end-2,end-1).*100,'-o')

g = 1.35;
H_PM = (0.21.*((UUvec).^2))./g; % Pierson-Moskowitz Sig Height https://wikiwaves.org/Ocean-Wave_Spectra

plot(UUvec(1:end-1),H_PM(1:end-1).*100,'--k')


legend('Wind Time Step (T) = 500','Wind Time Step (T) = 10000','Pierson-Moskowitz H_{sig}','location','best')