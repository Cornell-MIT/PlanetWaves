clc
clear
close all


% Table S1
Earth.rho_s = 2.65*1000;
Earth.rho = 1.0*1000;
Earth.nu = 1e-6;
Earth.g = 9.81;

Mars.rho_s = 2.9*1000;
Mars.rho = 1.0*1000;
Mars.nu = 1e-6;
Mars.g = 3.71;

Titan_warm.rho_s = 0.95*1000;
Titan_warm.rho = 0.67*1000;
Titan_warm.nu = 3e-7;
Titan_warm.g = 1.35;

Titan_cold.rho_s = 0.95*1000;
Titan_cold.rho = 0.54*1000;
Titan_cold.nu = 6e-7;
Titan_cold.g = 1.35;

% Table S3
Mars_DistalFan.slope = 0.003;
Mars_DistalFan.width = 27;

Mars_CentralFan.slope = 0.01;
Mars_CentralFan.width = 27;

[~,bedload_dominated] = calc_riverine_flux(Mars,Mars_DistalFan);
display_results(ab,bedload_dominated,'MARS','DISTAL FAN')
[~,bedload_dominated] = calc_riverine_flux(Mars,Mars_CentralFan);
display_results(bc,bedload_dominated,'MARS','CENTRAL FAN')


% Table S4
Jezero.slope = 0.03;
Jezero.width = 45;

% Table S5
Titan_saraswati_min.slope = 4e-4;
Titan_saraswati_min.width = 175;
Titan_saraswati_max.slope = 0.002;
Titan_saraswati_max.width = 700;

% Table S6
Titan_VidFlum_min.slope = 0.0011;
Titan_VidFlum_min.width = 100;
Titan_VidFlum_max.slope = 0.0015;
Titan_VidFlum_max.width = 175;





function display_results(suspended_load,bedload,planet,rivername)

        
    
    fprintf('For %s on %s:\n',rivername,planet)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('BEDLOAD DOMINATED')
    fprintf('D50 = %.1f cm\n',bedload.D50*100)
    fprintf('Q = %.0f m3/s\n',bedload.Q)
    fprintf('Qs = %.1e m3/s\n',bedload.Qs)
    fprintf('H = %.1f m\n',bedload.H)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('SUSPENDED LOAD DOMINATED')
    fprintf('D50 = %.1f cm\n',suspended_load.D50*100)
    fprintf('Q = %.0f m3/s\n',suspended_load.Q)
    fprintf('Qs = %.1e m3/s\n',suspended_load.Qs)
    fprintf('H = %.1f m\n',suspended_load.H)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')



end