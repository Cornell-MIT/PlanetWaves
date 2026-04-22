clc
clear
close all

% COMPARE WITH AYAKO'S RESULTS
addpath(genpath(fullfile('..','subroutines')));

Titan_warm.rho_s = 0.95*1000;
Titan_warm.rho = 0.67*1000;
Titan_warm.nu = 3e-7;
Titan_warm.g = 1.35;

Titan_cold.rho_s = 0.95*1000;
Titan_cold.rho = 0.54*1000;
Titan_cold.nu = 6e-7;
Titan_cold.g = 1.35;



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

width_range = [175, 700];
slope_range = [4e-4, 0.002];

figure;
hold on;

D50 = NaN(numel(width_range), numel(slope_range));

for j = 1:numel(slope_range)
    River.slope = slope_range(j);

    for i = 1:numel(width_range)
        River.width = width_range(i);

        [sand, gravel] = calc_riverine_flux(Titan_warm, River);

        D50(i,j) = sand.D50*100;
        Q(i,j)   = sand.Q;
        Qs(i,j)  = sand.Qs;
        H(i,j)   = sand.H;
    end

    % Plot D50 vs width for each slope as seperate line
    plot(width_range, H(:,j), '-o', 'DisplayName', sprintf('Slope = %.4g', River.slope));
end

yline(0.26,'--k','DisplayName','Sam Min')
yline(18,'--k','DisplayName','Sam Max')
%ylim([0 10])
xlabel('River Width');
ylabel('H');
legend('show','Location','best');
grid on;



% [susload,bedload_dominated] = calc_riverine_flux(Titan_cold,Titan_VidFlum_min);
% display_results(susload,bedload_dominated,'Titan (cold)','Vid Flumina (min)')
% [susload,bedload_dominated] = calc_riverine_flux(Titan_cold,Titan_VidFlum_max);
% display_results(susload,bedload_dominated,'Titan (cold)','Vid Flumina (max)')
% [susload,bedload_dominated] = calc_riverine_flux(Titan_warm,Titan_VidFlum_min);
% display_results(susload,bedload_dominated,'Titan (warm)','Vid Flumina (min)')
% [susload,bedload_dominated] = calc_riverine_flux(Titan_warm,Titan_VidFlum_max);
% display_results(susload,bedload_dominated,'Titan (warm)','Vid Flumina (max)')
% 
% 
% [susload,bedload_dominated] = calc_riverine_flux(Titan_cold,Titan_saraswati_min);
% display_results(susload,bedload_dominated,'Titan (cold)','Saraswati (min)')
% [susload,bedload_dominated] = calc_riverine_flux(Titan_cold,Titan_saraswati_max);
% display_results(susload,bedload_dominated,'Titan (cold)','Saraswati (max)')
% [susload,bedload_dominated] = calc_riverine_flux(Titan_warm,Titan_saraswati_min);
% display_results(susload,bedload_dominated,'Titan (warm)','Saraswati (min)')
% [susload,bedload_dominated] = calc_riverine_flux(Titan_warm,Titan_saraswati_max);
% display_results(susload,bedload_dominated,'Titan (warm)','Saraswati (max)')


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