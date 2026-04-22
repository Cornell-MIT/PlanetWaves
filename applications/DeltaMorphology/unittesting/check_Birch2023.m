clc
clear
close all

addpath(genpath(fullfile('..','subroutines')));

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

% [~,bedload_dominated] = calc_riverine_flux(Mars,Mars_DistalFan);
% display_results(NaN,bedload_dominated,'MARS','DISTAL FAN')
% [~,bedload_dominated] = calc_riverine_flux(Mars,Mars_CentralFan);
% display_results(NaN,bedload_dominated,'MARS','CENTRAL FAN')

% Table S4
Jezero.slope = 0.03;
Jezero.width = 45;

% [susload,bedload_dominated] = calc_riverine_flux(Mars,Jezero);
% display_results(susload,bedload_dominated,'Mars','JEZERO')


% Table S5
Titan_saraswati_min.slope = 4e-4;
Titan_saraswati_min.width = 175;
Titan_saraswati_max.slope = 0.002;
Titan_saraswati_max.width = 700;

% [susload,bedload_dominated] = calc_riverine_flux(Titan_cold,Titan_saraswati_min);
% display_results(susload,bedload_dominated,'Titan (cold)','Saraswati (min)')
% [susload,bedload_dominated] = calc_riverine_flux(Titan_cold,Titan_saraswati_max);
% display_results(susload,bedload_dominated,'Titan (cold)','Saraswati (max)')
% [susload,bedload_dominated] = calc_riverine_flux(Titan_warm,Titan_saraswati_min);
% display_results(susload,bedload_dominated,'Titan (warm)','Saraswati (min)')
% [susload,bedload_dominated] = calc_riverine_flux(Titan_warm,Titan_saraswati_max);
% display_results(susload,bedload_dominated,'Titan (warm)','Saraswati (max)')

% Table S6
Titan_VidFlum_min.slope = 0.0011;
Titan_VidFlum_min.width = 100;
Titan_VidFlum_max.slope = 0.0015;
Titan_VidFlum_max.width = 175;


[susload,bedload_dominated] = calc_riverine_flux(Titan_cold,Titan_VidFlum_min);
display_results(susload,bedload_dominated,'Titan (cold)','Vid Flumina (min)')
[susload,bedload_dominated] = calc_riverine_flux(Titan_cold,Titan_VidFlum_max);
display_results(susload,bedload_dominated,'Titan (cold)','Vid Flumina (max)')
[susload,bedload_dominated] = calc_riverine_flux(Titan_warm,Titan_VidFlum_min);
display_results(susload,bedload_dominated,'Titan (warm)','Vid Flumina (min)')
[susload,bedload_dominated] = calc_riverine_flux(Titan_warm,Titan_VidFlum_max);
display_results(susload,bedload_dominated,'Titan (warm)','Vid Flumina (max)')


%%%%
slopes = [0.0004,0.0013,0.002]; 
widths = [175:25:700];
j = 0;
for S = 1:numel(slopes) % slope
  j = j + 1;
    i = 0;
    for B = 1:numel(widths) %channel width
        i = i + 1;
        ThisRiver.slope = slopes(S);
        ThisRiver.width = widths(B);
        [susload,bedload] = calc_riverine_flux(Titan_warm,ThisRiver);
        myD50(j,i) = bedload.D50;
        myQ(j,i) = bedload.Q;
        myQs(j,i) = bedload.Qs;
        myD50_sand(j,i) = susload.D50;
        myQ_sand(j,i) = susload.Q;
        myQs_sand(j,i) = susload.Qs;
        myRe_p(j,i) = susload.Re_p;

    end
end

load("Ayako_SaraswatiFlumen_gravel.mat","D50_gravel","Q_gravel","Qs_gravel")

make_comparison_plot(myD50.*100,D50_gravel.*100,'the dimensional bed grain size (cm)','Saraswati Flumen, cold gravel')
make_comparison_plot(myQ,Q_gravel,'flow discharge','Saraswati Flumen, cold gravel')
make_comparison_plot(myQs,Qs_gravel,'sed flux','Saraswati Flumen, cold gravel')

load('Ayako_SaraswatiFlumen_sand.mat','D50_sand','Q_sand','Qs_sand','Rep')

make_comparison_plot(myD50_sand.*100,D50_sand.*100,'the dimensional bed grain size (cm)','Saraswati Flumen, cold sand')
make_comparison_plot(myQ_sand,Q_sand,'flow discharge','Saraswati Flumen, cold sand')
make_comparison_plot(myRe_p,Rep,'Particle Reynolds Number','Saraswati Flumen, cold sand')
make_comparison_plot(myQs_sand,Qs_sand,'sed flux','Saraswati Flumen, cold sand')




function make_comparison_plot(my_Vals,Ayako_Vals,name_val,regime)

    widths = [175:25:700];
    slopes = [0.0004,0.0013,0.002]; 
    clrs = winter(numel(slopes));

   
    figure;
    subplot(1,2,1)
    for i = 1:numel(slopes)
        plot(widths,Ayako_Vals(i,:),'--','LineWidth',3,'Color',clrs(i,:),'DisplayName',sprintf('S = %f, Ayako',slopes(i)))
        hold on
        plot(widths,my_Vals(i,:),'-','LineWidth',2,'Color',clrs(i,:),'DisplayName',sprintf('S = %f, Me',slopes(i)))
    end
    xlabel('channel width (m)')
    ylabel(name_val)
    title(regime)
    legend('show','Location','best')
    subplot(1,2,2)
    for i = 1:numel(slopes)
        plot(widths,round(Ayako_Vals(i,:)-my_Vals(i,:),2),'-','LineWidth',3,'Color',clrs(i,:))
        hold on
    end
    xlabel('channel width (m)')
    ylabel(name_val)
    if Ayako_Vals(i,:)-my_Vals(i,:) < 1e-2
        title(['no difference for ',regime])
    else
        title(['Differences of ',regime])
    end
    
    fprintf('[%s] Range of %s: %.2f - %.2f\n',regime,name_val,min(my_Vals,[],'all'),max(my_Vals,[],'all'))
end

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