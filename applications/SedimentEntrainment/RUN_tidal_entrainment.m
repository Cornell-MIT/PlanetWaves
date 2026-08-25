clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will produce the Shields diagram for tidal flow predictions
% within Titan's nearshore and straits within its lakes. These results are
% based on the tidal flow models of Tokano+2010, Tokano+2014, Vincent+2016,
% and Vincent+2018. It will also print to the console key values of interest
% for minimum grain sizes to mobilize as bedload and suspended load
% across a range of grain densities 
% Author: Una Schneck (ugschneck@gmail.com)
% Last Modified: 8/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TITAN GRAINS IN NEARSHORE AND STRAITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d50_titan = 6.35e-5:1e-6:0.1;                                                    % range for my grains [very fines gravel]
g = 1.352;                                                                 % titan gravity
rho_s = [800 940 1500];                                                    % grain density [organic-ice ice organic]

% PREVIOUS TIDAL MODEL PARAMETERS AND RESULTS
rho_strait = 550;                                                          % liquid density of straits kg/m3 (Vincent 2018)
rho_lake = 662;                                                            % liquid density of nearshore kg/m3 (Vincent 2016)
dyn_vis_lake = 1736e-6;                                                    % liquid viscosity of lake Pa.S (Lorenz 2010)
kin_vis_lake = dyn_vis_lake/rho_lake;                                      % kinematic viscosity of nearshore m2/s
kin_vis_strait = 3e-7;                                                     % kinematic ciscotiy strait m2/s (Vincent 2018)

fluid_strait = [rho_strait kin_vis_strait];
fluid_lake = [rho_lake kin_vis_lake];

max_depth_strait = 15;                                                     % m (V-shaped basin) (Vincent 2018)
depth_strait = max_depth_strait/2;                                         % use the average? the reported velocity would be a depth-averaged one
depth_lake = 3;                                                            % ~3 m (Vincent 2016, emailed to ask for better precision but not response yet

man_coef_max_strait = 0.03;                                                % Manning coefficient strait max Vincent 2018
man_coef_min_strait = 0.06;                                                % Manning coefficient strait min Vincent 2018
man_coef_lake = 0.03;                                                      % Manning coefficient nearshore Vincent 2016 

u_max_strait = 0.64;                                                       % Maximum flow velocity in strait m/s (Vincent 2018) table 3
u_min_strait = 0.12;                                                       % Minimum flow velocity in strait  m/s (Vincent 2018) table 3
u_max_lake = 0.046;                                                        % Maximum flow velocity in nearshore m/s (Vincent 2016)  
u_min_lake = 0.02;                                                         % Minimum flow velocity in lake m/s (Vincent 2016)

% possible liquid compositions across Titan's lakes
rho_liq_pos = [510.65 654.69];
nu_liq_pos = [3.1827e-7 1.6501e-6];

% CACULATE SHIELDS FOR RANGE IN NEARSHORE AND STRAITS
for i = 1:length(rho_s)

    % STRAITS
    Re_strait_max(i,:) = calc_particle_Reynolds(rho_s(i),rho_strait,g,d50_titan,kin_vis_strait);
    shields_strait_max(i,:) = calc_shields_number(rho_s(i),rho_strait,depth_strait,u_max_strait,d50_titan,man_coef_max_strait,g);

    Re_strait_min(i,:) = calc_particle_Reynolds(rho_s(i),rho_strait,g,d50_titan,kin_vis_strait);
    shields_strait_min(i,:) = calc_shields_number(rho_s(i),rho_strait,max_depth_strait,u_min_strait,d50_titan,man_coef_min_strait,g);

    % LAKES (NEARSHORE)
    % Re_lake_max(i,:) = calc_particle_Reynolds(rho_s(i),rho_lake, g, d50, kin_vis_lake);
    Re_lake_max(i,:) = calc_particle_Reynolds(rho_s(i),rho_lake,g,d50_titan,kin_vis_lake);
    shields_lake_max(i,:) = calc_shields_number(rho_s(i), rho_lake, depth_lake, u_max_lake, d50_titan, man_coef_lake,g);

    %Re_lake_min(i,:) = calc_particle_Reynolds(rho_s(i),rho_lake, g, d50, kin_vis_lake);
    Re_lake_min(i,:) = calc_particle_Reynolds(rho_s(i),rho_lake,g,d50_titan,kin_vis_lake);
   
    shields_lake_min(i,:) =  calc_shields_number(rho_s(i), rho_lake, depth_lake, u_min_lake, d50_titan, man_coef_lake,g);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE ENTRAINMENT AND SUSPENSION THRESHOLDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Earth values 
e_rho = 998.57;                                                            % water at 18C
e_rho_s = 2620;                                                            % quartz sand
e_g = 9.81;                                                                % earth gravity
e_nu = 1.0533e-6;                                                          % water viscosity (m2/s) water at 18C
d50_earth = 6.35e-6:1e-6:0.2;                                                    % range of median grain size for thresholds [very fines gravel] (making it longer and more continuous than grains)

% particle entrainment thresholds
Re_p = calc_particle_Reynolds(e_rho_s,e_rho,e_g,d50_earth,e_nu);
threshold_entrainment1 = calculate_Parker_shields_curve(Re_p);
threshold_entrainment2 = calc_original_shields_curve(Re_p);

% particle suspension thresholds
[suspension,s2,s3,s4] = calc_suspension_threshold(e_rho_s,e_rho,e_g,d50_earth,e_nu);

% making the plot
idx = unique(min(round(logspace(0, log10(length(Re_p)), 500)), length(Re_p)));
x = Re_p(idx).^2;
y1 = threshold_entrainment1(idx);
y2 = threshold_entrainment2(idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pretty colors for plot
HEXCOLORS = {'#fc8d59','#91bfdb','#FADA5E'};
for i = 1:numel(HEXCOLORS)
    str = HEXCOLORS{i};
    mycolors(i,:) = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
end

figure('Name','Tidal Flow Competence');
grid on;
hold on;
set(gca,'YScale','log');
set(gca,'XScale','log')
xlabel('Dimensionless Particle Size, $D_*$','Interpreter','latex')
ylabel('Dimensionless Shear Stress, $\theta$','Interpreter','latex')
xlim([10^-2.5 10^11]);
ylim([10^-4 10^3]);

% (1) plot the entrainment thresholds
fill([x, fliplr(x)], [y1, fliplr(y2)], [0.7 0.7 0.7], 'EdgeColor', 'none','HandleVisibility', 'off'); 

% (2) plot the suspension thresholds
mask_valid = suspension > threshold_entrainment1;
Re_p_valid = Re_p(mask_valid);
susp_valid = suspension(mask_valid);
s2_valid = s2(mask_valid);
if length(Re_p_valid) < 2
    warning('Not enough points for find threshold');
    return
end
% find where threshold for entrainment and suspension meet so can plot only
% the part of the suspension threshold that is above the entrainment 
% threshold
x1 = Re_p(suspension > threshold_entrainment1).^2;
y1 = suspension(suspension > threshold_entrainment1);
x2 = Re_p(s2 > threshold_entrainment1).^2;
y2 = s2(s2 > threshold_entrainment1);
[x_common, ia, ib] = intersect(x1, x2);                                    % need common x-values for fill
if isempty(x_common)
    x_fill = linspace(max(min(x1), min(x2)), min(max(x1), max(x2)), 200);
    y1i = interp1(x1, y1, x_fill, 'linear', 'extrap');
    y2i = interp1(x2, y2, x_fill, 'linear', 'extrap');
else
    x_fill = x_common;
    y1i = y1(ia);
    y2i = y2(ib);
end
idx = round(logspace(log10(2), log10(length(x_fill)), 500));
idx = unique(min(idx, length(x_fill)));
idx = idx(:); 
x_fill_sus = [x_fill(idx), fliplr(x_fill(idx))];
y_fill_sus = [y1i(idx), fliplr(y2i(idx))];
fill(x_fill_sus, y_fill_sus,[0.7 0.7 0.7], 'EdgeColor', 'none','HandleVisibility', 'off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPORT KEY VALUES TO CONSOLE (and add vertical lines to Shields diagram)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grain_labels = {'Light Organics Grains [$\rho$ = 800 kg/m$^3$]', 'Water Ice Grains [$\rho$ = 940 kg/m$^3$]','Dense Organics Grains [$\rho$ = 1500 kg/m$^3$]'};
grain_names = {'Light organics', 'Water ice', 'Dense organics'};            % without rho


% get the index for the dimensionless grain size
get_log_indices = @(N) unique(min(round(logspace(0, log10(N), 500)), N));
for i = 1:numel(rho_s)
    % draw vertical lines corresponding to reported values on Shields diagram to check
    idx = get_log_indices(length(Re_strait_max(i,:)));
    x = Re_strait_max(i, idx).^2;
    y1 = shields_strait_max(i, idx);
    y2 = shields_strait_min(i, idx);
    x2 = [x, fliplr(x)];
    yfill = [y1, fliplr(y2)];
    fill(x2, yfill, 'r', 'FaceAlpha', 0.5, 'LineStyle', 'none', 'FaceColor', mycolors(i,:), 'HandleVisibility','off');
end


for i = 1:numel(rho_s)
    idx = get_log_indices(length(Re_lake_max(i,:)));
    x = Re_lake_max(i, idx).^2;
    y1 = shields_lake_max(i, idx);
    y2 = shields_lake_min(i, idx);
    x2 = [x, fliplr(x)];
    yfill = [y1, fliplr(y2)];

    fill(x2, yfill, 'k', 'FaceAlpha', 0.5, 'LineStyle', 'none', 'FaceColor', mycolors(i,:),'DisplayName', grain_labels{i});
end

legend('show', 'Interpreter', 'latex');

% report effect of liquid on shields parameter
Re_liq_min = calc_particle_Reynolds(rho_s(2),rho_liq_pos(1),g,d50_titan,nu_liq_pos(1));
Re_liq_max = calc_particle_Reynolds(rho_s(2),rho_liq_pos(2),g,d50_titan,nu_liq_pos(2));
shields_liq_min = calc_shields_number(rho_s(2),rho_liq_pos(1),depth_strait,u_max_strait,d50_titan,man_coef_max_strait,g);
shields_liq_max = calc_shields_number(rho_s(2),rho_liq_pos(2),depth_strait,u_max_strait,d50_titan,man_coef_max_strait,g);
plot(Re_liq_min.^2,shields_liq_min,'--b','DisplayName','Minimum Density Liquid')
plot(Re_liq_max.^2,shields_liq_max,'--r','DisplayName','Maximum Density Liquid')
x1 = Re_liq_min.^2;
x2 = Re_liq_max.^2;
y1 = shields_liq_min;
y2 = shields_liq_max;
% Find the overlap
x_min = max(min(x1), min(x2));
x_max = min(max(x1), max(x2));
x_common = linspace(x_min, x_max, 100); 
y1_interp = interp1(x1, y1, x_common, 'linear');
y2_interp = interp1(x2, y2, x_common, 'linear');
avg_diff = mean(y2_interp - y1_interp);
fprintf('Average difference (y2 - y1) over overlapping x: %.4f\n', avg_diff);

% report key D50 values
fprintf('\n\n');
fprintf('====================================================================\n');
fprintf(' MAXIMUM TRANSPORTABLE GRAIN SIZES\n');
fprintf('====================================================================\n');

% STRAIT
fprintf('\n');
fprintf('######################## STRAIT ##############################\n');
fprintf('Maximum velocity = %.3f m/s\n',u_max_strait);
fprintf('Depth used       = %.3f m\n',depth_strait);
fprintf('Fluid density    = %.1f kg/m^3\n',rho_strait);
fprintf('Kinematic visc.  = %.3e m^2/s\n',kin_vis_strait);
fprintf('################################################################\n');

for i = 1:length(rho_s)

    % Critical entrainment threshold
    theta_crit = calculate_Parker_shields_curve(Re_strait_max(i,:));
    % Suspension threshold
    [theta_susp,~,~,~] = calc_suspension_threshold(rho_s(i),rho_strait,g,d50_titan,kin_vis_strait);
    % Maximum-flow Shields stress
    theta_flow = shields_strait_max(i,:);
    % BEDLOAD
    mobile = theta_flow >= theta_crit;
    if any(mobile)
        D_bedload = max(d50_titan(mobile));
    else
        D_bedload = NaN;
    end
    % SUSPENSION
    suspended = (theta_flow >= theta_susp) & mobile;
    if any(suspended)
        D_suspended = max(d50_titan(suspended));
    else
        D_suspended = NaN;
    end
   
    fprintf('\n%s (rho = %.0f kg/m^3)\n', grain_names{i},rho_s(i));

    if isnan(D_bedload)
        fprintf('\tMaximum bedload grain:     NONE\n');
    else
        fprintf('\tMaximum bedload grain:     %.6f m  = %.3f mm\n', D_bedload,D_bedload*1000);
    end

    if isnan(D_suspended)
        fprintf('\tMaximum suspended grain:   NONE\n');
    else
        fprintf('\tMaximum suspended grain:   %.6f m  = %.3f mm\n', D_suspended,D_suspended*1000);
    end

    if ~isnan(D_bedload) && ~isnan(D_suspended)
        fprintf('\tBedload-only range:        %.3f--%.3f mm\n', D_suspended*1000,D_bedload*1000);
    end

    strait_results(i).D_bedload = D_bedload;
    strait_results(i).D_suspended = D_suspended;

end

% NEARSHORE
fprintf('\n');
fprintf('######################## LAKE ###############################\n');
fprintf('Maximum velocity = %.3f m/s\n',u_max_lake);
fprintf('Depth             = %.3f m\n',depth_lake);
fprintf('Fluid density     = %.1f kg/m^3\n',rho_lake);
fprintf('Kinematic visc.   = %.3e m^2/s\n',kin_vis_lake);
fprintf('################################################################\n');

for i = 1:length(rho_s)
    % particle reynolds number
    Re_lake = calc_particle_Reynolds(rho_s(i),rho_lake,g,d50_titan,kin_vis_lake);
    % Critical entrainment threshold
    theta_crit = calculate_Parker_shields_curve(Re_lake);
    % Suspension threshold 
    [theta_susp,~,~,~] = calc_suspension_threshold(rho_s(i),rho_lake,g,d50_titan,kin_vis_lake);
    % Maximum-flow lake Shields stress
    theta_flow = shields_lake_max(i,:);
    % bedload indices
    mobile = theta_flow >= theta_crit;
    if any(mobile)
        D_bedload = max(d50_titan(mobile));
    else
        D_bedload = NaN;
    end
    % suspended load indices
    suspended = (theta_flow >= theta_susp) & mobile;

    if any(suspended)
        D_suspended = max(d50_titan(suspended));
    else
        D_suspended = NaN;
    end
    fprintf('\n%s (rho = %.0f kg/m^3)\n',grain_names{i},rho_s(i));

    if isnan(D_bedload)
        fprintf('\tMaximum bedload grain:     NONE\n');
    else
        fprintf('\tMaximum bedload grain:     %.6f m  = %.3f mm\n',D_bedload,D_bedload*1000);
    end

    if isnan(D_suspended)
        fprintf('\tMaximum suspended grain:   NONE\n');
    else
        fprintf('\tMaximum suspended grain:   %.6f m  = %.3f mm\n', D_suspended,D_suspended*1000);
    end

    if ~isnan(D_bedload) && ~isnan(D_suspended)
        fprintf('\tBedload-only range:        %.3f--%.3f mm\n', D_suspended*1000,D_bedload*1000);
    end

    lake_results(i).D_bedload = D_bedload;
    lake_results(i).D_suspended = D_suspended;
    lake_results(i).Re = Re_lake;

end


% DRAW LINES ON SHIELDS DIAGRAM (STRAIT)
for i = 1:length(rho_s)

    D = strait_results(i).D_bedload;
    if ~isnan(D)
        Re_D = calc_particle_Reynolds(rho_s(i),rho_strait,g,D,kin_vis_strait);

        x_D = Re_D^2;

        xline(x_D,'--','Color',mycolors(i,:), 'LineWidth',2, 'HandleVisibility','off');
        text(x_D, 2e-4, sprintf('%.2f mm',D*1000), 'Color',mycolors(i,:), 'FontSize',10, 'FontWeight','bold', 'Rotation',90,  'HorizontalAlignment','left', 'VerticalAlignment','bottom', 'Clipping','on');
    end

    D = strait_results(i).D_suspended;
    if ~isnan(D)

        Re_D = calc_particle_Reynolds(rho_s(i),rho_strait,g,D,kin_vis_strait);

        x_D = Re_D^2;

        xline(x_D,':', 'Color',mycolors(i,:), 'LineWidth',2,'HandleVisibility','off');
        text(x_D, 2e-4, sprintf('%.2f mm',D*1000), ...
            'Color',mycolors(i,:), 'FontSize',10, 'FontWeight','bold', 'Rotation',90, 'HorizontalAlignment','left', 'VerticalAlignment','bottom', 'Clipping','on');
    end
end


% DRAW LINES ON SHIELDS DIAGRAM (NEARSHORE)
for i = 1:length(rho_s)
    D = lake_results(i).D_bedload;
    if ~isnan(D)
        Re_D = calc_particle_Reynolds(rho_s(i),rho_lake,g,D,kin_vis_lake);
        x_D = Re_D^2;
        xline(x_D,'--', 'Color',mycolors(i,:),  'LineWidth',1.5, 'HandleVisibility','off');
    end

    D = lake_results(i).D_suspended;
    if ~isnan(D)
        Re_D = calc_particle_Reynolds(rho_s(i),rho_lake,g,D,kin_vis_lake);
        x_D = Re_D^2;
        xline(x_D,':', ...
            'Color',mycolors(i,:), 'LineWidth',1.5, 'HandleVisibility','off');
    end
end

% FIND GRAIN SIZE WHERE ENTRAINMENT AND SUSPENSION THRESHOLDS INTERSECT
Re_p = calc_particle_Reynolds(e_rho_s,e_rho,e_g,d50_earth,e_nu);

theta_entrain = calculate_Parker_shields_curve(Re_p);
[theta_susp,~,~,~] = calc_suspension_threshold(e_rho_s,e_rho,e_g,d50_earth,e_nu);

% Difference between the two thresholds
diff_theta = theta_susp - theta_entrain;

% Find sign changes
idx_cross = find(diff_theta(1:end-1).*diff_theta(2:end) <= 0);

for j = 1:length(idx_cross)

    k = idx_cross(j);

    % Interpolate in log(D50) for better accuracy at crossover
    logD_cross = interp1(diff_theta(k:k+1), log10(d50_earth(k:k+1)), 0);

    D_cross = 10^logD_cross;

    Re_cross = calc_particle_Reynolds(e_rho_s,e_rho,e_g,D_cross,e_nu);

    theta_cross = calculate_Parker_shields_curve(Re_cross);

    fprintf('\nThreshold intersection:\n');
    fprintf('\tD50 = %.6e m = %.4f mm = %.2f microns\n', D_cross,D_cross*1000,D_cross*1e6);
    fprintf('\tRe_p = %.4f\n',Re_cross);
    fprintf('\tD_*^2 = %.4f\n',Re_cross^2);
    fprintf('\tShields = %.5f\n',theta_cross);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- internal functions -------------------------------------------- %
function [s1,s2,s3,s4] = calc_suspension_threshold(rho_s,rho,g,D50,nu)

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
    u_star = settling_velocity;%;./(P*kappa);
    suspension_threshold = (rho.*(u_star.^2))./((rho_s - rho)*g.*D50);

    % figure;
    % plot(D_star,W_star,'--')
    % hold on;
    % 
    % 
    % set(gca,'YScale','log');
    % set(gca,'XScale','log')


    bagnold_sus = 0.64.*((rho.*(settling_velocity.^2))./((rho_s - rho).*g.*D50)); % bagnold1966
    parker_sus = (rho.*(settling_velocity.^2))./((rho_s - rho)*g.*D50); % from online slides of parker (https://studylib.net/doc/9626909/1d-sediment-transport-morphodynamics-with)
    for i = 1:numel(D_star)
        REP = sqrt(D_star);
        if REP(i) > 1 && REP(i) < 32
            K = 4*(REP(i)^(-2/3));
        elseif REP(i) > 32
            K = 0.4;
        else
            K = NaN;
        end
        van_rijn_1984(i) = sqrt(K)*parker_sus(i);
    end
    % plot(sqrt(D_star),bagnold_sus,'--r')
    % plot(sqrt(D_star),parker_sus,'--g')
    % plot(sqrt(D_star),van_rijn_1984,'--b')
    % legend('Dietrich','Bagnold','Parker','Van Rijn')
    % 
    % xlabel('$D_*$','Interpreter','latex')
    % ylabel('$\tau$','Interpreter','latex')

    s1 = suspension_threshold;
    s2 = bagnold_sus;
    s3 = parker_sus;
    s4 = van_rijn_1984;

    % title('suspension thresholds')
end

function tau_c_star = calc_original_shields_curve(s_star)
    tau_c_star = 0.105.*(s_star).^(-0.3) + 0.045.* exp(-35.* (s_star).^(-0.59));
end

function critical_shields = calculate_Parker_shields_curve(Re_particle)
    critical_shields = 0.5.*(0.22.*(Re_particle.^(-0.6)) + 0.06.*(10.^(-7.7.*(Re_particle.^(-0.6)))));
end

function Re_p = calc_particle_Reynolds(rho_s,rho,g,D50,nu)
    spec_weight = (rho_s/rho) - 1;
    Re_p = sqrt(spec_weight*g.*(D50.^3))./nu;
end

function D50 = calc_D50_from_Re_p(rho_s, rho, g, nu, Re_p)
    spec_weight = (rho_s / rho) - 1;
    D50 = ((Re_p .* nu).^2 ./ (spec_weight * g)).^(1/3);
end

function theta = calc_shields_number(rho_s,rho,depth,u,D50,Manning,g)
    %theta = (rho*Manning*(abs(u)^2))./((rho_s - rho)*g.*D50.*(depth^(1/3)));
    kappa = 0.4;
    z0 = D50/30;
    u_star = (kappa*abs(u))./log((0.37*depth)./(z0));
    theta = (rho*(u_star.^2))./((rho_s - rho)*g.*D50);
end

