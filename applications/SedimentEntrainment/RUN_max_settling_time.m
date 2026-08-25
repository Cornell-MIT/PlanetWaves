clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will calculate the maximum settling time for grains
% i.e., how long does the turbidity signal persist if all other motion
% ceases after the wave is passed
% This script will produce a plot of settling time (in hours) versus
% median grain size for different grain density and liquid composition
% Author: Una Schneck (ugschneck@gmail.com)
% Last Modified: 8/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath(fullfile('..','..','planetwaves')))

lakes = {'Titan-CH3H8N2','Titan-OntarioLacus','Titan-LigeiaMare','Titan-CH4N2'};
lakecolors = {'#F2D1C9','#E086D3','#8332AC','#462749'};
icecolor = {'#42F2F7','#46ACC2','#498C8A','#4B6858'};                      % ice colors (blues)
organicicecolor = {'#ffa5ab','#da627d','#a53860','#450920'};               % organic-ice (reds)
organiccolor =  {'#f8ed62','#e9d700','#dab600','#a98600'};                 % organics (yellows)

rho_s = [800 940 1500];                                                    % grain densities kg/m3 (light organic; ice; dense organic)
d50 = 6.35e-5:1e-5:0.1;                                                    % median grain sizes m [fines gravel]

d = [1 7.1 10];                                                            % depth of interest (7.1 corresponds to depth for average summer storm conditions)

for dd = 2%1:numel(d)
    depth = d(dd);
    figure('Name','Settling times for grains through liquid column');

for s = 1:numel(rho_s)
    if s == 1
        mycolor = organicicecolor;
    elseif s == 2
        mycolor = icecolor;
    elseif s == 3
        mycolor = organiccolor;
    end

    for c = 1:numel(lakes)
        
        
        [Planet, ~, ~, ~, ~] = initalize_model(lakes{c}, 0, 0, 0, [0 0]);

        settling_time = depth./calc_settling_velocity(rho_s(s),Planet.rho_liquid,Planet.gravity,d50,Planet.nu_liquid);
        plot(d50,settling_time./(60*60),'LineWidth',3,'Color',mycolor{c},'DisplayName',[lakes{c} ' (rho = ' num2str(rho_s(s)) ')'])
        hold on
        set(gca,'YScale','log')
        set(gca,'XScale','log')
    end
end
grid on
xlabel('D$_{50}$ (m)','Interpreter','latex','FontSize',20)
ylabel('Settling time (hr)','Interpreter','latex','FontSize',20)
grid on;
box on;
title(sprintf('Settling through column of %0.1f m',depth),'Interpreter','latex','FontSize',20)
ax = gca;
ax.FontSize = 18;
%legend('show')
end

function settling_velocity = calc_settling_velocity(rho_s,rho,g,D50,nu)

    P = 0.8;%1.2;
    kappa = 0.4;
    R = (rho_s/rho) - 1;
    D_star = (R*g.*(D50.^3))./(nu^2);
    log_D_star = log10(D_star);  
    log_W_star = -3.76715 + ...
                 1.92944 .* log_D_star - ...
                 0.09815 .* (log_D_star.^2) - ...
                 0.00575 .* (log_D_star.^3) + ...
                 0.00056 .* (log_D_star.^4);
    W_star = 10.^(log_W_star);

    settling_velocity = ((W_star.*R.*g.*nu)).^(1/3);

end