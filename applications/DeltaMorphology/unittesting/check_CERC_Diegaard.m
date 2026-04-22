clc
clear
close all


H_wave = [0.8 0.8 0.8];
L_wave = [25 25 25];
T_wave = [4 4 4];

planet = 'Earth';

if strcmp(planet,'Titan')
    g = 1.352;
    rho_s = 940;
    rho = 588.15; 
    nu = 8.084e-7;
elseif strcmp(planet,'Earth')
    g = 9.81;
    rho_s = 2650;
    rho = 1000;
    nu = 1.2e-6;
end

alpha0 = wrapToPi(deg2rad(0:359));

beta = 1e-3;
D50 = 6.35e-5; % d50 = 6.35e-5:1e-5:0.1;  % [fines gravel]

Qs_Diegaard = calc_Diegaard_wave_flux(H_wave, L_wave, rho_s, rho, g, D50, nu, beta, alpha0);
Qs_CERC = CERC(H_wave,T_wave,rho,rho_s,g,rad2deg(alpha0));

[Qs_max,i_max] = max(Qs_Diegaard);
[Qs_min,i_min] = min(Qs_Diegaard);
[QsC_max,ic_max] = max(Qs_CERC);
[QsC_min,ic_min] = min(Qs_CERC);

figure;
plot(rad2deg(alpha0),Qs_Diegaard./max(Qs_Diegaard),'ob','DisplayName','Diegaard')
hold on;
plot(rad2deg(alpha0),Qs_CERC./max(Qs_CERC),'.r','DisplayName','CERC')
grid on;
xline(rad2deg(alpha0(i_max)),'-b',sprintf('max angle: %.0f',rad2deg(alpha0(i_max))),'LabelVerticalAlignment','middle','HandleVisibility','off')
xline(rad2deg(alpha0(i_min)),'-b',sprintf('min angle: %.0f',rad2deg(alpha0(i_min))),'LabelVerticalAlignment','middle','HandleVisibility','off')
xline(rad2deg(alpha0(ic_max)),'-r',sprintf('max angle: %.0f',rad2deg(alpha0(ic_max))),'LabelVerticalAlignment','middle','HandleVisibility','off')
xline(rad2deg(alpha0(ic_min)),'-r',sprintf('min angle: %.0f',rad2deg(alpha0(ic_min))),'LabelVerticalAlignment','middle','HandleVisibility','off')

xlabel('relative approach angle [deg]')
ylabel('Qs/Qsmax')
legend('show','Location','best')


figure;
plot(rad2deg(alpha0),Qs_Diegaard,'ob','DisplayName','Diegaard')
hold on;
plot(rad2deg(alpha0),Qs_CERC,'.r','DisplayName','CERC')
grid on;
xlabel('relative approach angle [deg]')
ylabel('Qs')
legend('show','Location','best')

