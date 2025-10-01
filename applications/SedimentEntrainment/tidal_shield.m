clc
clear
close all

% PLOT OF TIDAL FLOW ENTRAINMENT FOR TITAN STRAITS AND NEARSHORE FLOWS

rho_liq_pos = [510.65 654.69];
nu_liq_pos = [3.1827e-7 1.6501e-6];
% --- Constants ---- %
% Earth values 
e_rho = 998.57; % water at 18C
e_rho_s = 2620; % quartz 
e_g = 9.81; % earth gravity
e_nu = 1.0533e-6; % m2/s water at 18C

% --- Plotting Elements ---- %

HEXCOLORS = {'#fc8d59','#91bfdb','#FADA5E'};
for i = 1:numel(HEXCOLORS)
    str = HEXCOLORS{i};
    mycolors(i,:) = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
end

d50 = [6.35e-6:1e-6:0.2];  % [very fines gravel]


% --- Particle Entrainment ---- %
Re_p = calc_particle_Reynolds(e_rho_s,e_rho,e_g,d50,e_nu);
threshold_entrainment1 = calculate_Parker_shields_curve(Re_p);
threshold_entrainment2 = calc_original_shields_curve(Re_p);


% --- Particle Suspension --- %
[suspension,s2,s3,s4] = calc_suspension_threshold(e_rho_s,e_rho,e_g,d50,e_nu);


% --- Plotting entrainment threshold ---- %
nPoints = 500;
N = length(Re_p);
idx = unique(min(round(logspace(0, log10(N), nPoints)), N));
x = Re_p(idx).^2;
y1 = threshold_entrainment1(idx);
y2 = threshold_entrainment2(idx);

figure;

% plot(x, y1, '-k', 'LineWidth', 1, 'HandleVisibility', 'off');
% plot(x, y2, '-k', 'HandleVisibility', 'off');
fill([x, fliplr(x)], [y1, fliplr(y2)], [0.7 0.7 0.7], 'EdgeColor', 'none','HandleVisibility', 'off'); % entrainment threshold

grid on;
hold on;



% Plotting suspension threshold  --- %
mask_valid = suspension > threshold_entrainment1;
Re_p_valid = Re_p(mask_valid);
susp_valid = suspension(mask_valid);
s2_valid = s2(mask_valid);


N = length(Re_p_valid);
if N < 2
    warning('Too few valid points after thresholding.');
    return
end

idx = round(logspace(log10(2), log10(N), nPoints));
idx = unique(min(idx, N));
idx = idx(:); 


x = Re_p_valid(idx).^2;
y1 = susp_valid(idx);
y2 = s2_valid(idx);
x = x(:); y1 = y1(:); y2 = y2(:);
y1 = max(y1, threshold_entrainment1);
y2 = max(y2, threshold_entrainment1);

x1 = Re_p(suspension > threshold_entrainment1).^2;
y1 = suspension(suspension > threshold_entrainment1);

x2 = Re_p(s2 > threshold_entrainment1).^2;
y2 = s2(s2 > threshold_entrainment1);

[x_common, ia, ib] = intersect(x1, x2); % need common x-values for fill

if isempty(x_common)
    x_fill = linspace(max(min(x1), min(x2)), min(max(x1), max(x2)), 200);
    y1i = interp1(x1, y1, x_fill, 'linear', 'extrap');
    y2i = interp1(x2, y2, x_fill, 'linear', 'extrap');
else
    x_fill = x_common;
    y1i = y1(ia);
    y2i = y2(ib);
end


nPoints = 500;
N = length(x_fill);
idx = round(logspace(log10(2), log10(N), nPoints));
idx = unique(min(idx, N));
idx = idx(:); 

x_fill_sus = [x_fill(idx), fliplr(x_fill(idx))];
y_fill_sus = [y1i(idx), fliplr(y2i(idx))];



fill(x_fill_sus, y_fill_sus,[0.7 0.7 0.7], 'EdgeColor', 'none','HandleVisibility', 'off'); % suspension threshold

set(gca,'YScale','log');
set(gca,'XScale','log')
axis([10^-1 10^7 10^-3 10^3])
xlabel('Dimensionless Particle Size, $D_*$','Interpreter','latex')
ylabel('Dimensionless Shear Stress, $\theta$','Interpreter','latex')


% ----- Plot tidal flow shields parameters ------ %

d50 = [6.35e-5:1e-6:0.1];  % [very fines gravel]

g = 1.352;
rho_s = [800 940 1500]; % [organic-ice ice organic]

% PREVIOUS TIDAL MODEL PARAMETERS AND RESULTS
rho_strait = 550; % kg/m3 (Vincent 2018)
rho_lake = 662; % kg/m3 (Vincent 2016)

dyn_vis_lake = 1736e-6; % Pa.S (Lorenz 2010)
kin_vis_lake = dyn_vis_lake/rho_lake; % m2/s
kin_vis_strait = 3e-7; % m2/s (Vincent 2018)

fluid_strait = [rho_strait kin_vis_strait];
fluid_lake = [rho_lake kin_vis_lake];

max_depth_strait = 15; % m (V-shaped basin) (Vincent 2018)
depth_strait = max_depth_strait/2; % use the average? the reported velocity would be a depth-averaged one
depth_lake = 3; % ~3 m (Vincent 2016, emailed to ask for better precision but not response yet

man_coef_max_strait = 0.03; % Vincent 2018
man_coef_min_strait = 0.06; % Vincent 2018
man_coef_lake = 0.03; % Vincent 2016 

u_max_strait = 0.64; % m/s (Vincent 2018) table 3
u_min_strait = 0.12; % m/s (Vincent 2018) table 3
u_max_lake = 0.046; % m/s (Vincent 2016)  
u_min_lake = 0.02; % m/s (Vincent 2016)


for i = 1:length(rho_s)

    % STRAITS
    Re_strait_max(i,:) = calc_particle_Reynolds(rho_s(i),rho_strait,g,d50,kin_vis_strait);
    shields_strait_max(i,:) = calc_shields_number(rho_s(i),rho_strait,depth_strait,u_max_strait,d50,man_coef_max_strait,g);

    Re_strait_min(i,:) = calc_particle_Reynolds(rho_s(i),rho_strait,g,d50,kin_vis_strait);
    shields_strait_min(i,:) = calc_shields_number(rho_s(i),rho_strait,max_depth_strait,u_min_strait,d50,man_coef_min_strait,g);

    % LAKES (NEARSHORE)
    Re_lake_max(i,:) = calc_particle_Reynolds(rho_s(i),rho_lake, g, d50, kin_vis_lake);
    shields_lake_max(i,:) = calc_shields_number(rho_s(i), rho_lake, depth_lake, u_max_lake, d50, man_coef_lake,g);


    Re_lake_min(i,:) = calc_particle_Reynolds(rho_s(i),rho_lake, g, d50, kin_vis_lake);
    shields_lake_min(i,:) =  calc_shields_number(rho_s(i), rho_lake, depth_lake, u_min_lake, d50, man_coef_lake,g);

end

nPoints = 500;
get_log_indices = @(N) unique(min(round(logspace(0, log10(N), nPoints)), N));
for i = 1:3
    N = length(Re_strait_max(i,:));
    idx = get_log_indices(N);

    x = Re_strait_max(i, idx).^2;
    y1 = shields_strait_max(i, idx);
    y2 = shields_strait_min(i, idx);

    x2 = [x, fliplr(x)];
    yfill = [y1, fliplr(y2)];

    fill(x2, yfill, 'r', 'FaceAlpha', 0.5, 'LineStyle', 'none', 'FaceColor', mycolors(i,:), 'HandleVisibility','off');
end

labels = { ...
    'Light Organics Grains [$\rho$ = 800 kg/m$^3$]', ...
    'Water Ice Grains [$\rho$ = 940 kg/m$^3$]', ...
    'Dense Organics Grains [$\rho$ = 1500 kg/m$^3$]'};

for i = 1:3
    N = length(Re_lake_max(i,:));
    idx = get_log_indices(N);

    x = Re_lake_max(i, idx).^2;
    y1 = shields_lake_max(i, idx);
    y2 = shields_lake_min(i, idx);

    x2 = [x, fliplr(x)];
    yfill = [y1, fliplr(y2)];

    fill(x2, yfill, 'k', 'FaceAlpha', 0.5, 'LineStyle', 'none', 'FaceColor', mycolors(i,:),'DisplayName', labels{i});
end

legend('show', 'Interpreter', 'latex');
xlim([10^-2.5 10^11]);
ylim([10^-4 10^3]);

Re_liq_min = calc_particle_Reynolds(rho_s(2),rho_liq_pos(1),g,d50,nu_liq_pos(1));
Re_liq_max = calc_particle_Reynolds(rho_s(2),rho_liq_pos(2),g,d50,nu_liq_pos(2));

shields_liq_min = calc_shields_number(rho_s(2),rho_liq_pos(1),depth_strait,u_max_strait,d50,man_coef_max_strait,g);
shields_liq_max = calc_shields_number(rho_s(2),rho_liq_pos(2),depth_strait,u_max_strait,d50,man_coef_max_strait,g);

plot(Re_liq_min.^2,shields_liq_min,'--b','DisplayName','Minimum Density Liquid')
plot(Re_liq_max.^2,shields_liq_max,'--r','DisplayName','Maximum Density Liquid')


x1 = Re_liq_min.^2;
x2 = Re_liq_max.^2;
y1 = shields_liq_min;
y2 = shields_liq_max;

% 1. Find the overlap
x_min = max(min(x1), min(x2));
x_max = min(max(x1), max(x2));

% 2. Create a common x-vector in the overlap region
x_common = linspace(x_min, x_max, 100); % finer resolution if needed

% 3. Interpolate both datasets to x_common
y1_interp = interp1(x1, y1, x_common, 'linear');
y2_interp = interp1(x2, y2, x_common, 'linear');

% 4. Compute average difference
avg_diff = mean(y2_interp - y1_interp);

% Display result
fprintf('Average difference (y2 - y1) over overlapping x: %.4f\n', avg_diff);
% Re_lake_max_print = calc_particle_Reynolds(rho_s(3),rho_lake, g, d50, kin_vis_lake);
% 
% Re = calc_particle_Reynolds(rho_s(3),rho_lake, g,  1.2533e-04, kin_vis_lake);
% D50 = calc_D50_from_Re_p(rho_s(3), rho_lake, g, kin_vis_lake, 0.7);


%D50 = calc_D50_from_Re_p(rho_s(3), rho_strait, g, kin_vis_strait, sqrt(6.6438))
%D50 = calc_D50_from_Re_p(rho_s(3), rho_strait, g, kin_vis_strait, sqrt(971.144))
%D50 = calc_D50_from_Re_p(rho_s(2), rho_strait, g, kin_vis_strait, sqrt(303.014))
%D50 = calc_D50_from_Re_p(rho_s(2), rho_strait, g, kin_vis_strait, sqrt(10156.4))
%D50 = calc_D50_from_Re_p(rho_s(1), rho_strait, g, kin_vis_strait, sqrt(2951.05))
% D50 = calc_D50_from_Re_p(rho_s(1), rho_strait, g, kin_vis_strait, sqrt(28394))
D50 = calc_D50_from_Re_p(rho_s(1), rho_strait, g, kin_vis_strait, sqrt(0.93))

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
    plot(sqrt(D_star),bagnold_sus,'--r')
    plot(sqrt(D_star),parker_sus,'--g')
    plot(sqrt(D_star),van_rijn_1984,'--b')
    legend('Dietrich','Bagnold','Parker','Van Rijn')

    xlabel('$D_*$','Interpreter','latex')
    ylabel('$\tau$','Interpreter','latex')

    s1 = suspension_threshold;
    s2 = bagnold_sus;
    s3 = parker_sus;
    s4 = van_rijn_1984;

    title('suspension thresholds')
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

% 
% % % plot(Re_strait_max(1,:),shields_strait_max(1,:),'Color',mycolors(1,:),'LineWidth',5,'DisplayName','Ice-Organic Grains [$\rho$ = 800 kg/m$^3$]');    % MAX STRAIT FOR ICE-ORGANICS
% % % plot(Re_strait_max(2,:),shields_strait_max(2,:),'Color',mycolors(2,:),'LineWidth',5,'DisplayName','Ice Grains [$\rho$ = 940 kg/m$^3$]');            % MAX STRAIT FOR ICE
% % % plot(Re_strait_max(3,:),shields_strait_max(3,:),'Color',mycolors(3,:),'LineWidth',5,'DisplayName','Organic Grains [$\rho$ = 1500 kg/m$^3$]');        % MAX STRAIT FOR ORGANICS
% % % 
% % % 
% % % plot(Re_strait_min(1,:),shields_strait_min(1,:),'Color',mycolors(1,:),'LineWidth',5,'HandleVisibility','off');              % MIN STRAIT FOR ICE-ORGANICS
% % % plot(Re_strait_min(2,:),shields_strait_min(2,:),'Color',mycolors(2,:),'LineWidth',5,'HandleVisibility','off');              % MIN STRAIT FOR ICE
% % % plot(Re_strait_min(3,:),shields_strait_min(3,:),'Color',mycolors(3,:),'LineWidth',5,'HandleVisibility','off');              % MIN STRAIT FOR ORGANICS
% % 
% % % FILL IN ERROR BAR FOR STRAIT
% % curve1 = shields_strait_max(1,:); % MAX STRAIT ICE-ORGANIC
% % curve2 = shields_strait_min(1,:); % MIN STRAIT ICE-ORGANIC
% % x2 = [Re_strait_max(1,:),fliplr(Re_strait_max(1,:))];
% % inbetween = [curve1,fliplr(curve2)];
% % f2 = fill(x2,inbetween,'r','FaceAlpha',1,'LineStyle','none','FaceColor',mycolors(1,:),'HandleVisibility','off');
% % 
% % curve1 = shields_strait_max(2,:); % MAX STRAIT ICE 
% % curve2 = shields_strait_min(2,:); % MIN STRAIT ICE
% % x2 = [Re_strait_max(2,:),fliplr(Re_strait_max(2,:))];
% % inbetween = [curve1,fliplr(curve2)];
% % f1 = fill(x2,inbetween,'r','FaceAlpha',1,'LineStyle','none','FaceColor',mycolors(2,:),'HandleVisibility','off');
% % 
% % 
% % curve1 = shields_strait_max(3,:); % MAX STRAIT ORGANIC
% % curve2 = shields_strait_min(3,:); % MIN STRAIT ORGANIC
% % x2 = [Re_strait_max(3,:),fliplr(Re_strait_max(3,:))];
% % inbetween = [curve1,fliplr(curve2)];
% % f3 = fill(x2,inbetween,'r','FaceAlpha',1,'LineStyle','none','FaceColor',mycolors(3,:),'HandleVisibility','off');
% % 
% % % plot(Re_lake_max(1,:),shields_lake_max(1,:),'-','Color',mycolors(1,:),'LineWidth',5,'HandleVisibility','off'); % MAX STRAIT FOR ICE-ORGANICS
% % % plot(Re_lake_max(2,:),shields_lake_max(2,:),'-','Color',mycolors(2,:),'LineWidth',5,'HandleVisibility','off'); % MAX STRAIT FOR ICE
% % % plot(Re_lake_max(3,:),shields_lake_max(3,:),'-','Color',mycolors(3,:),'LineWidth',5,'HandleVisibility','off'); % MAX STRAIT FOR ICE
% % % 
% % % 
% % % plot(Re_lake_min(1,:),shields_lake_min(1,:),'-','Color',mycolors(1,:),'LineWidth',5,'HandleVisibility','off'); % MIN STRAIT FOR ICE-ORGANICS
% % % plot(Re_lake_min(2,:),shields_lake_min(2,:),'-','Color',mycolors(2,:),'LineWidth',5,'HandleVisibility','off'); % MIN STRAIT FOR ICE
% % % plot(Re_lake_min(3,:),shields_lake_min(3,:),'-','Color',mycolors(3,:),'LineWidth',5,'HandleVisibility','off'); % MIN STRAIT FOR ICE
% % 
% % % FILL IN ERROR BAR FOR STRAIT
% % curve1 = shields_lake_max(1,:); % MAX LAKE ICE-ORGANIC
% % curve2 = shields_lake_min(1,:); % MIN LAKE ICE-ORGANIC
% % x2 = [Re_lake_max(1,:),fliplr(Re_lake_max(1,:))];
% % inbetween = [curve1,fliplr(curve2)];
% % f2 = fill(x2,inbetween,'r','FaceAlpha',1,'LineStyle','none','FaceColor',mycolors(1,:),'DisplayName','Ice-Organic Grains [$\rho$ = 800 kg/m$^3$]');
% % 
% % curve1 = shields_lake_max(2,:); % MAX STRAIT ICE 
% % curve2 = shields_lake_min(2,:); % MIN STRAIT ICE
% % x2 = [Re_lake_max(2,:),fliplr(Re_lake_max(2,:))];
% % inbetween = [curve1,fliplr(curve2)];
% % f1 = fill(x2,inbetween,'r','FaceAlpha',1,'LineStyle','none','FaceColor',mycolors(2,:),'DisplayName','Ice Grains [$\rho$ = 940 kg/m$^3$]');
% % 
% % 
% % curve1 = shields_lake_max(3,:); % MAX LAKE ORGANIC
% % curve2 = shields_lake_min(3,:); % MIN LAKE ORGANIC
% % x2 = [Re_lake_max(3,:),fliplr(Re_lake_max(3,:))];
% % inbetween = [curve1,fliplr(curve2)];
% % f3 = fill(x2,inbetween,'r','FaceAlpha',1,'LineStyle','none','FaceColor',mycolors(3,:),'DisplayName','Organic Grains [$\rho$ = 1500 kg/m$^3$]'); 
% % 
% % legend('show','Interpreter','latex')
% % xlim([10^-2 10^7])
% % ylim([10^-3 10^3])
% 
% %exportgraphics(gcf, 'Tidal_Shield.pdf', 'ContentType', 'vector');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
