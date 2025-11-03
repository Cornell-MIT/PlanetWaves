clc
clear
close all

% Deriving a non-dimensional growth curve akin to Kahma+1981 which found 
% an eqn relating the non-dim energy to the non-dim fetch.

% final result: [u,x] -> [Hs,Tp]

% addpath(fullfile('..','..','planetwaves'))  
% addpath(fullfile('..','..','planetwaves','pre_analysis'))
% 
% 
% % MODEL INPUTS
test_speeds = [2:10];
grid_resolution = [0.1*1000 0.1*1000];
% downstream_num = 5;
% crossstream_num = 20;
% 
% time_to_run = 60*10;  
% wind_direction = 0;  
% zDep = 100.*ones(crossstream_num,downstream_num);
% buoy_loc = [round(downstream_num/2) round(crossstream_num/2)];    
% planet_to_run = 'Earth';
% [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
% Model.gridX = grid_resolution(1);                                              
% Model.gridY = grid_resolution(2);  
% 
% %make_input_map(Planet,Model,Wind)
% 
% for i = 1:numel(test_speeds)
% 
%     Wind.speed = test_speeds(i);
%     Model = calc_cutoff_freq(Planet,Model,Wind);
% 
%     [myHsig{i}, htgrid{i}, wn_e_spectrum, ~ , ~ , ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
%     if ~isempty(wn_e_spectrum{end})
%         energy{i} = squeeze(sum(wn_e_spectrum{end}.E(Model.long,Model.lat,:,:),4));
%         wn{i} = squeeze(sum(wn_e_spectrum{end}.k(Model.long,Model.lat,:,:),4));
%         cg{i} = squeeze(sum(wn_e_spectrum{end}.cg(Model.long,Model.lat,:,:),4));
%     end
% end
% save('Earth_Fetch3.mat','myHsig','htgrid','wn','energy','cg')
% make_plots(Planet,Model,Wind,test_speeds,myHsig, htgrid,energy,wn)

load('Earth_Fetch.mat','htgrid')

clr = winter(numel(test_speeds));
figure
tiledlayout("horizontal")
for jj = 1:numel(htgrid{1})
    clf
    nexttile;
for u = 1:numel(htgrid)
    if ~isempty(htgrid{u}{jj})
        crossstream(u,:) = htgrid{u}{jj}(3,:);
        downstream(u,:) = htgrid{u}{jj}(:,10);

        downstream(u,1:2) = NaN; 
        downstream(u,end-1:end) = NaN;
        crossstream(u,1:2) = NaN; 
        crossstream(u,end-1:end) = NaN;

        fetch = grid_resolution(1).*(1:numel(downstream(u,:)));

        plot(fetch./1000,downstream(u,:),'Color',clr(u,:),'DisplayName',num2str(test_speeds(u)))
        hold on
    end
end
xlabel('fetch km')
ylabel('wave height m')
ylim([0 1])

g = 9.81;
wind_speeds = test_speeds;

X_nd = zeros(numel(wind_speeds), numel(fetch));    
E_nd = zeros(numel(wind_speeds), numel(fetch));   

for i = 1:numel(wind_speeds)
    for j = 1:numel(fetch)

        u = wind_speeds(i);
        x = fetch(j); % Non-dimensional fetch
        Hs = downstream(i,j); % Non-dimensional energy

        X_nd(i,j) = calc_nondim_fetch(x, u , g);

        E_nd(i,j) = calc_nondim_energy(Hs, u, g);
    end
end

nexttile
for u = 1:numel(test_speeds)
    plot(X_nd(u,:), E_nd(u,:), 'o','Color',clr(u,:),'MarkerFaceColor',clr(u,:),'DisplayName',num2str(test_speeds(u)))
    hold on
end
xlabel('Non-dimensional Fetch')
ylabel('Non-dimensional Energy')
title(['Timestep: ',num2str(jj)])
grid on

X_THEORY = logspace(2,5);
E_JONSWAP = (1.6e-7).*X_THEORY;
E_BOTHNIAN = (3.6e-7).*X_THEORY;
E_PM_MAX = 3.64e-3.*ones(size(X_THEORY));
E_CERC = (5e-3).*((tanh(0.0125.*(X_THEORY.^0.42))).^2);
E_LAKEONTARIO = (8.415e-7).*(X_THEORY.^0.76);
E_DOBSON = (12.7e-7).*(X_THEORY.^0.75);

plot(X_THEORY,E_JONSWAP,'--k','DisplayName','JONSWAP, fetch limited')
plot(X_THEORY,E_BOTHNIAN,'--k','DisplayName','Bothnian Sea, fetch limited')
plot(X_THEORY,E_CERC,'--k','DisplayName','CERC, fetch limit')
plot(X_THEORY,E_LAKEONTARIO,'--k','DisplayName','Lake Ontario, fetch limit')
plot(X_THEORY,E_DOBSON,'--k','DisplayName','Dobson, fetch limit')
plot(X_THEORY,E_PM_MAX,'--k','DisplayName','PM, fully developed')

set(gca,'XScale','log')
set(gca,'YScale','log')
xlim([10^2 10^5])
ylim([10^-5 10^-1])

% if jj == 1
%     gif('nondim_thru_time.gif');
% else
%     gif;
% end
end

%legend('show')
function X = calc_nondim_fetch(x,u,g)
    % calculates non-dimensional fetch for a constant wind speed
    % INPUT
    %   x = fetch [1 x num_fetch]
    %   u = wind speed [1 x 1]
    %   g = gravity [1 x 1]
    % OUTPUT
    %   X = non-dimensional fetch [1 x num_fetch]

    X = (g.*x)./(u^2);

end

function e = calc_nondim_energy(Hs,u,g)
    % calculates non dimensional energy for a constant fetch and changing wind speed
    % INPUT
    %   Hs = significant wave height [1 x num_wind_speeds]
    %   u  = wind speed [1 x num_wind_speeds]
    %   g  = gravity [1 x 1]
    % OUTPUT
    %   e  = non-dimensional energy [1 x num_wind_speeds]
   
    h = Hs./4;
    e = ((h.^2).*(g^2))./(u.^4);

end