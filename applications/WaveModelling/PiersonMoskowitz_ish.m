% clc
% clear
% close all
% 
% % compute the Pierson-Moskowitz curve for wave height
% 
% addpath(genpath(fullfile('..','..','planetwaves')))
% 
% save_file = ['PM_curve_fits_', datestr(datetime("today")), '.mat']; % name of file with saved data
% 
% % Model
% planet_to_run  = 'Earth';
% test_speeds    = 10;            % m/s
% time_to_run    = 60*10;         % 10 hours
% wind_direction = 0;
% buoy_loc       = [5 5];
% grid_resolution = [10*1000 1*1000];
% 
% zDep = 100 .* ones(10,10);
% 
% % Initialize model
% [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run, time_to_run, wind_direction, zDep, buoy_loc);
% Model.gridX = grid_resolution(1);
% Model.gridY = grid_resolution(2);
% 
% % sweep across ratios
% rat = [0.25 0.5 1 2 4];
% 
% % parameters to sweep across
% paramList = { ...
%     'rho_liquid', ...
%     'nu_liquid', ...
%     'nua', ...
%     'gravity', ...
%     'surface_temp', ...
%     'surface_press', ...
%     'surface_tension', ...
%     'kgmolwt', ...
%     };
% 
% 
% % SWEEP AWAY, BABY!!! (all but wind speed and fetch)
% results = struct();
% for k = 1:numel(paramList)
% 
%     paramName = paramList{k};
% 
%     fprintf('\n*******Running power-law sweep for %s*************\n', paramName);
% 
%     [n,rat,Hsig] = powerlaw_wave_sweep(paramName, rat,Planet, Model, Wind, Uniflow, Etc, test_speeds);
% 
%     results.(paramName).n    = n;
%     results.(paramName).rat  = rat;
%     results.(paramName).Hsig = Hsig;
% 
% end
% 
% % % sweep for wind speed
% Urat = rat;
% Wind.speed = test_speeds;
% fprintf('\n************* Running power-law sweep for wind speed *************\n');
% [n_U, Urat, HU] = powerlaw_wind_sweep(Urat, Planet, Model, Wind, Uniflow, Etc);
% fprintf('Hs \\alpha U^{%.3f}\n', n_U);

% sweep for fetch
fprintf('\n ************* Running power-law sweep for fetch *************\n');
[n_fetch, fetch_m, Hfetch] = powerlaw_fetch_sweep(Planet, Model, Wind, Uniflow, Etc, test_speeds, 5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FINAL SUMMARY OF SCALINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n============================================\n');
fprintf('  FINAL POWER-LAW SUMMARY (Hs scaling)\n');
fprintf('============================================\n');

% parameters besides wind and fetch
for k = 1:numel(paramList)
    pname = paramList{k};
    fprintf('%-20s : %8.4f\n', pname, results.(pname).n);
end

% Wind speed
fprintf('%-20s : %8.4f\n', 'wind_speed', n_U);
% Fetch
fprintf('%-20s : %8.4f\n', 'fetch', n_fetch);
fprintf('============================================\n');


%save(save_file,"results")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POWER LAW SWEEPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n, rat, Hsig] = powerlaw_wave_sweep(paramName, rat, Planet, Model, Wind, Uniflow, Etc, test_speeds)
%   Sweep over a planet parameter and extracts power-law exponent
%       Hs ∝ (paramName)^n

    if ~isfield(Planet, paramName)
        error('Planet has no field named "%s"', paramName);
    end
    
    % hold reference value
    p0 = Planet.(paramName);
    pa0 = Planet.rhoa;
    % initialize
    Hsig = zeros(size(rat));
    
    % SWEEP
    for i = 1:numel(rat)
    
        % Reset parameter 
        Planet.(paramName) = p0 * rat(i);

        if strcmp(paramName,'surface_press') || strcmp(paramName, 'kgmolwt')
            Planet.rhoa =  pa0* rat(i);
        elseif strcmp(paramName,'surface_temp')
           Planet.rhoa = pa0 / rat(i);
        end

        % Reset wind
        Wind.speed = test_speeds;
        % Recompute cutoff freq 
        Model = calc_cutoff_freq(Planet, Model, Wind);
        % Run wave model
        [myHsig,~,~,~,~,~,~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);
        % store final wave height
        Hsig(i) = myHsig(end);
    
    end
    
    % normalize by typical Earth conditions (ratio = 1)
    idx = (rat == 1);
    if ~any(idx)
        error('rat must include 1 for normalization');
    end
    
    Href = Hsig(idx);
    X = rat;
    Y = Hsig / Href;
    
    % fit to power law
    logX = log(X);
    logY = log(Y);
    p = polyfit(logX, logY, 1);
    n = p(1);
    % Plot fit
    figure; hold on; grid on;
    plot(X, Y, 'ko', 'MarkerFaceColor','k')
    plot(X, exp(polyval(p, logX)), 'r--', 'LineWidth',1.5)
    set(gca,'XScale','log','YScale','log')
    xlabel(['$', paramName, '/', paramName, '_0$'], 'Interpreter','latex')
    ylabel('$H_s / H_{s,0}$', 'Interpreter','latex')

    title_name = replace(paramName, "_", "\_");
    title(['Power law: $H_s \propto ', title_name, '^{', num2str(n,'%.f'), '}$'], 'Interpreter','latex')
    drawnow;
end


function [n, Urat, Hsig] = powerlaw_wind_sweep( Urat, Planet, Model, Wind, Uniflow, Etc)
%   Sweeps wind speed and extracts power-law exponent
%
%   Hs ∝ U_10^n

    % hold reference wind speed
    U0 = Wind.speed;
    % initalize 
    Hsig = zeros(size(Urat));
    
    % SWEEP
    for i = 1:numel(Urat)
    
        Wind.speed = U0 * Urat(i);
    
        % recompute cutoff freq as needed
        Model = calc_cutoff_freq(Planet, Model, Wind);
        % Run wave model
        [myHsig,~,~,~,~,~,~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);
        % final wave height
        Hsig(i) = myHsig(end);
    
    end
    
    % normalize to a typical wind speed (e.g., 10 m/s)
    idx = (Urat == 1);
    if ~any(idx)
        error('Urat must include 1 for normalization');
    end
    Href = Hsig(idx);
    X = Urat;
    Y = Hsig / Href;
    
    % fit power law
    logX = log(X);
    logY = log(Y);
    p = polyfit(logX, logY, 1);
    n = p(1);
    % power law
    figure; hold on; grid on;
    plot(X, Y, 'ko', 'MarkerFaceColor','k')
    plot(X, exp(polyval(p, logX)), 'b--', 'LineWidth',1.5)
    set(gca,'XScale','log','YScale','log')
    xlabel('$U / U_0$', 'Interpreter','latex')
    ylabel('$H_s / H_{s,0}$', 'Interpreter','latex')
    title(['Wind scaling: $H_s \propto U^{',num2str(n,'%.f'), '}$'], 'Interpreter','latex')

end


function [n, fetch_m, H_vs_fetch] = powerlaw_fetch_sweep(Planet, Model, Wind, Uniflow, Etc, test_speed, col_idx)
%   Sweep over fetches to extract power law 
%       Hs ∝ X^n
  
    if nargin < 7
        col_idx = round(Model.LatDim/2);   % default cross-fetch column (center of basin)
    end
    
    % model inputs
    Wind.speed = test_speed;
    Model = calc_cutoff_freq(Planet, Model, Wind);
    make_input_map(Planet, Model, Wind);
    % run model
    [~, htgrid, ~, ~, ~, ~, ~] = makeWaves(Planet, Model, Wind, Uniflow, Etc);
    % Extract final height grid
    plot_grid = htgrid{end}';
    plot_grid(isnan(plot_grid)) = 0;
    % Fetch axis (meters)
    n_rows = size(plot_grid,1);
    fetch_m = ((1:n_rows) * Model.gridY)';
    
    % H vs fetch at fixed cross-fetch location
    H_vs_fetch = plot_grid(:, col_idx);
    % clean up
    valid_idx = 1:round(Model.LonDim/2);
    fetch_m = fetch_m(valid_idx);
    H_vs_fetch = H_vs_fetch(valid_idx);
    
    % fit to power law
    logF = log(fetch_m);
    logH = log(H_vs_fetch);
    coeffs = polyfit(logF, logH, 1);
    n = coeffs(1);
    a = exp(coeffs(2));
    % plot power law
    figure; hold on; grid on;
    loglog(fetch_m, H_vs_fetch, 'ko', 'MarkerFaceColor','k')
    loglog(fetch_m, a * fetch_m.^n, 'r--', 'LineWidth',1.5)
    xlabel('Fetch X [m]', 'Interpreter','latex')
    ylabel('$H_s$ [m]', 'Interpreter','latex')
    title(sprintf('Fetch scaling: $H_s \propto \\,X^{%.f}$', a, n),'Interpreter','latex')
    fprintf('Estimated fetch law: H \propto \\X^{%.f}\n', a, n);

end

