% clc
% clear
% close all
% 
% addpath(fullfile('..','..','planetwaves'))  
% addpath(fullfile('..','..','planetwaves','pre_analysis'))  
% 
% % MAKE PIERSON-MOSKOWITZ-ish CURVE FOR CHANGING PLANETARY CONDITIONS FOR HSIG AND UMIN
% 
% planet_to_run = 'Earth';
% test_speeds = 10;
% time_to_run = 60*10;  
% wind_direction = pi/2;  
% buoy_loc = [3 40];    
% grid_resolution = [1*1000 0.01*1000];
% zDep = 100.*ones(80,10);
% 
% [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
% 
% Model.gridX = grid_resolution(1);                                              
% Model.gridY = grid_resolution(2);   
% 
% 
% % (1)
% %Earth = Planet.kgmolwt;
% rat = [0.25 0.5 1 2 4];
% %var = Earth.*rat;
% %test_speeds = test_speeds.*rat;
% 
% 
% figure;
% time_evolve_ax = axes;
% grid on;
% legend('show', 'Location', 'northwest','interpreter','latex');
% title(['Waves on',' ',Planet.name,' at ',num2str(Planet.surface_press/1000), ' kPa'],'interpreter','latex');
% xlabel('model time step [$\Delta$ t]','interpreter','latex')
% ylabel('significant wave height [m]','interpreter','latex')
% hold on;
% 
% 
% for i = 1:numel(test_speeds)
% 
%                 Wind.speed = test_speeds(i);
% % (2)
%                 %Planet.kgmolwt = var(i);
%                 Model = calc_cutoff_freq(Planet,Model,Wind);
% 
%                 make_input_map(Planet,Model,Wind)
%                 [myHsig,htgrid, ~, ~ , ~ , ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
%                 plot(time_evolve_ax,1:numel(myHsig),myHsig,'-','DisplayName',num2str(Wind.speed))
%                 Hsig(i) = myHsig(end);
% 
% end
% 
% plot_grid = htgrid{end}';
% plot_grid(isnan(plot_grid)) = 0;
% plot_alpha_data = ones(size(plot_grid));
% plot_alpha_data(plot_grid==0) = 0;
% 
% figure;
% imagesc(plot_grid)
% hold on;
% xline(5, 'r--', 'LineWidth', 2);
% colorbar;
% set(gca,'YDir','reverse');      
% colormap jet;
% colorbar;
% 
% 
% new_xtick = get(gca, 'XTick') * Model.gridX / 1000;
% new_ytick = get(gca, 'YTick') * Model.gridY / 1000;
% set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%d', x), new_xtick, 'UniformOutput', false));
% set(gca, 'YTickLabel', arrayfun(@(y) sprintf('%d', y), new_ytick, 'UniformOutput', false));
% 
% 
% 
% 
% n_rows = size(plot_grid, 1);
% fetch_km = ((1:n_rows) * Model.gridY) / 1000;
% 
% H_vs_fetch = plot_grid(:, 5);
% 
% figure;
% plot(fetch_km,H_vs_fetch)
% xlabel('fetch km')
% ylabel('H m')
% 
% 
fetch_m = fetch_km * 1000;  
fetch_m = fetch_m';
valid_idx = fetch_m > 0 & H_vs_fetch > 0;
fetch_m = fetch_m(valid_idx);
H_vs_fetch = H_vs_fetch(valid_idx);

logF = log(fetch_m);
logH = log(H_vs_fetch);
coeffs = polyfit(logF, logH, 1);
b = coeffs(1);
log_a = coeffs(2);
a = exp(log_a);
fprintf('Estimated power law: H = %.4f * F^{%.4f}\n', a, b);

figure;
scatter(fetch_m, H_vs_fetch, 'b.');
hold on;
% % my_rat = Hsig./Hsig(rat == 1)
% % 
% % % fit power law
% % X = rat; 
% % Y = my_rat;
% logX = log(X);
% logY = log(Y);
% fitResult = polyfit(logX, logY, 1);
% n = fitResult(1);
% fprintf('The best fit exponent n is: %f\n', n);
% 





