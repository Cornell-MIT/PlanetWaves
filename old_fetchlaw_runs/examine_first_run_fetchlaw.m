clc
clear 
close all

load('T_500.mat','UUvec','ht') % ran for Time = 500

figure
plot(UUvec(2:end-1),ht(2:end-2,end-1).*100,'-o')
hold on
xlabel('wind speed [m/s]')
ylabel('sig H [cm]')
grid on

load('T_10000.mat','UUvec','ht') % ran for Time = 10000


plot(UUvec(2:end-1),ht(2:end-2,end-1).*100,'-o')

g = 1.35;
H_PM = (0.21.*((UUvec).^2))./g; % Pierson-Moskowitz Sig Height https://wikiwaves.org/Ocean-Wave_Spectra

plot(UUvec(2:end-1),H_PM(2:end-1).*100,'--k')

xlim([0 3.5])



addpath 'C:\Users\schne\OneDrive\Desktop\Main\Work\MIT\Github_Repos\ShorelineShaping\waves\user_defined_functions'

for i = 1:length(0.4:0.1:3.3)
    
    u1 = 0.4 + 0.1*i;
    u2 = u1;
    [~,H0,~,~] = make_wave_COE([u1 u2],1e6,100,1.2,465, 0.0031/1e4,'Titan',0);
    
    u(i) = u1;
    H(i) = H0(1)
end

plot(u,H.*100,'ok')

legend('Wind Time Step (T) = 500','Wind Time Step (T) = 10000','Pierson-Moskowitz H_{sig}','H Monochromatic Energy Conserve','location','best')
