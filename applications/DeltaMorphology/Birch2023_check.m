clc
clear
close all


% Table S1
Earth.rho_s = 2.65*1000; % kg/m3
Earth.rho = 1*1000; % kg/m3
Earth.nu = 0.01/10000; % m2/s
Earth.g = 9.81;

Mars.rho_s = 2.9*1000; 
Mars.rho = 1*1000;
Mars.nu = 0.01/10000;
Mars.g = 3.71;

Titan_cold.rho_s = 0.95*1000; 
Titan_cold.rho = 0.67*1000;
Titan_cold.nu = 0.003/10000;
Titan_cold.g = 1.35;

Titan_warm.rho_s = 0.95*1000;
Titan_warm.rho = 0.54*1000;
Titan_warm.nu = 0.006/10000;
Titan_warm.g = 1.35;

pos_planets = {Earth,Mars,Titan_cold,Titan_warm};
% labels = {'','','cold bedload','warm bedload'};
planet_label = {'Earth','Mars','cold Titan','warm Titan'};
river_label = {'','','Vid Flumina','Vid Flumina'};

saraswati.coldbedload_D50 = [1.6 61]./100;
saraswati.warmbedload_D50 = [1.5 56]/100;
saraswati.width = [175 700];
saraswati.slope = [4e-4 0.002];

vidflumina.coldbedload_D50 = [3.8 10]/100;
vidflumina.warmbedload_D50 = [3.5 9.4]/100;
vidflumina.susload_D50 = 6.35e-5;
vidflumina.width = [100 175];
vidflumina.slope = [0.0011 0.0015];

for p = 3:4% cold Titan
    bedload_Q_range = [];
    bedload_Qs_range = [];
    bedload_H_range = [];
    bedload_S_range = [];
    bedload_B_range = [];
    susload_Q_range = [];
    susload_Qs_range = [];
    susload_H_range = [];
    susload_S_range = [];
    susload_B_range = [];

    planet = pos_planets{p};


    if p == 3
        %D50 = saraswati.coldbedload_D50; 
        gravel_D50 = vidflumina.coldbedload_D50;
    end
    if p == 4
        %D50 = saraswati.warmbedload_D50;
        gravel_D50 = vidflumina.warmbedload_D50;
    end
    %river_width = saraswati.width; % Saraswati Flumen
    river_width = vidflumina.width;
    
    for d = 1:numel(gravel_D50)
        for b = 1:numel(river_width)
            [sand_river,gravel_river] = riverine_flux(planet.rho_s,planet.rho,planet.nu,planet.g,gravel_D50(d),vidflumina.susload_D50,river_width(b),NaN);


            bedload_Q_range = [bedload_Q_range gravel_river.Q];
            bedload_Qs_range = [bedload_Qs_range gravel_river.Qs];
            bedload_H_range = [bedload_H_range gravel_river.H];
            bedload_S_range = [bedload_S_range gravel_river.S];
            bedload_B_range = [bedload_B_range gravel_river.B];

            susload_Q_range = [susload_Q_range sand_river.Q];
            susload_Qs_range = [susload_Qs_range sand_river.Qs];
            susload_H_range = [susload_H_range sand_river.H];
            susload_S_range = [susload_S_range sand_river.S];
            susload_B_range = [susload_B_range sand_river.B];
        end
    end

    disp('------')
    disp('------')
    fprintf('%s at %s:\n',planet_label{p},river_label{p})
    disp('BEDLOAD')
    fprintf('Q: %g - %g  m3/s\n',min(bedload_Q_range),max(bedload_Q_range))
    fprintf('Q_s: %g - %g m3/s\n',min(bedload_Qs_range),max(bedload_Qs_range))
    fprintf('H: %g - %g m\n',min(bedload_H_range),max(bedload_H_range))
    fprintf('S: %g - %g\n',min(bedload_S_range),max(bedload_S_range))
    fprintf('B: %g - %g m\n',min(bedload_B_range),max(bedload_B_range))
    disp('------')
    disp('SUSPENDED LOAD')
    fprintf('Q: %g - %g m3/s\n',min(susload_Q_range),max(susload_Q_range))
    fprintf('Qs: %g - %g m3/s\n',min(susload_Qs_range),max(susload_Qs_range))
    fprintf('H: %g - %g m\n',min(susload_H_range),max(susload_H_range))
    fprintf('S: %g - %g \n',min(susload_S_range),max(susload_S_range))
    fprintf('B: %g - %g m\n',min(susload_B_range),max(susload_B_range))


    
end



