clc
clear
close all

% use latext font
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

% good font sizes for exporting
set(groot, 'defaultAxesLabelFontSizeMultiplier', 1);
set(groot, 'defaultAxesFontSize', 18);  % tic font size
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.11);
set(groot, 'defaultAxesFontName', 'default');  % default LaTeX font

label_font_size =  20;


addpath(fullfile('..','past_runs'))
addpath(fullfile('..','..','planetwaves'))
addpath(fullfile('..','..','planetwaves','post_analysis'))

addpath(fullfile('..','..','data/Titan/TAMwTopo/'))
addpath(fullfile('..','..','data/Titan/TitanLakes/Bathymetries/bathtub_bathy'))


lakecolors = gray(5);

% Results from Titan_Waves03_43
load('TAM_OL_winds','mag_wind','angle_wind');

u = round(mag_wind,1);
psi = round(angle_wind,1);

edges =  0:0.5:max(u) + 0.5; % slice every 0.5 m/s
[counts, ~] = histcounts(u, edges);  
frac_time = counts / numel(u);

load('../past_runs/Waves03_43.mat','test_speeds','time_to_run','wind_direction','zDep','buoy_loc','lakes','Planet','Model','H0','L0','T','C0','Cg0')


lake_slope = 2e-3;

min_depth = 10^-4;
alpha_0 = 0;
d50 = [6.35e-5:1e-5:0.1];  % [fines gravel]

d = max(max(Model.bathy_map)):-lake_slope:min_depth;
d50 = linspace(d50(1),d50(end),100);
rho_s = [800 940 1500]; % [organic-ice ice organic]
max_steepness = 1/7;

ol = 2; lm = 3;
org = 1; ice = 2; fluffy = 3;

max_d0_shoal = [];
% CALCULATE ENTRAINMENT DEPTH FOR SAND AND GRAVEL
for c = 1:numel(lakes)

    [Planet, Model, Wind, Uniflow, Etc] = initalize_model(lakes{c}, time_to_run, wind_direction, zDep, buoy_loc);

    %smooth artifacts in wave model (introduces max of ~ 0.1 error)
    %figure
    %plot(test_speeds,H0(1,:),'-','DisplayName','original H')
    %hold on
    H0(c,:) = smooth_waves(H0(c,:));
    %plot(test_speeds,H0(1,:),'--','DisplayName','smoothed H')
    %figure;
    %plot(test_speeds,T(1,:),'-','DisplayName','original T')
    %hold on
    T(c,:) = smooth_waves(T(c,:));
    %plot(test_speeds,T(1,:),'--','DisplayName','smoothed T')
    %figure;
    %plot(test_speeds,L0(1,:),'-','DisplayName','original L0')
    L0(c,:) = smooth_waves(L0(c,:));
    %hold on
    %plot(test_speeds,L0(1,:),'--','DisplayName','smoothed L0')


    for u = 1:numel(test_speeds)

            % shoaling waves
            L_shoal = L0(c,u) * sqrt(tanh(4*(pi^2)*d/(T(c,u)^2*Planet.gravity)));
            d_L = d ./ L_shoal;
            alpha = asin(tanh((2*pi.*d_L).*sin(alpha_0)));
            C_shoal = C0(c,u) * tanh(2*pi*d_L);
            n = 0.5 * (1 + ((4*pi*d_L)./sinh(4*pi*d_L)));
            Cg_shoal = n .* C_shoal;
            KR = sqrt(cos(alpha_0)./cos(alpha));
            KS = sqrt(Cg0(c,u)./Cg_shoal);
            H_shoal = KR .* KS .* H0(c,u);
            % orbital size, velocity, breaking wave condition
            d0_shoal = (H_shoal / 2) .* (1 ./ sinh(2*pi*d_L));
            if test_speeds(u) == 0.7 || test_speeds(u) == 0.8
                max_d0_shoal = [max_d0_shoal max(d0_shoal,[],'omitnan')];
            end
            um_shoal = (H_shoal / 2) .* ((Planet.gravity * T(c,u)) ./ L_shoal) .* (1 ./ cosh(2*pi*d_L));
            break_frac = H_shoal ./ L_shoal;
            break_frac(break_frac >= max_steepness) = NaN;

             % meshgrid to vectorizing for computational speed
            [S, A] = meshgrid(rho_s, d50);  % [S,A] : [numel(d50), numel(rho_s)]

            shields = zeros(length(d), size(S, 1), size(S, 2));  % size : [numel(d), numel(d50), numel(rho_s)]
            KM_crash = zeros(length(d), size(A, 1), size(A, 2));  % size : [numel(d), numel(d50), numel(rho_s)]

            for z = 1:length(d)
                shields(z, :, :) = (Planet.rho_liquid * (um_shoal(z)^2)) ./ ((S - Planet.rho_liquid) * Planet.gravity .* A);
                KM_crash(z, :, :) = 0.3 .* sqrt(d0_shoal(z) ./ A);  
            end

            for s = 1:numel(rho_s)
                for a = 1:numel(d50)
                    % across all depths                    
                    shield_across_depth = squeeze(shields(:, a, s));  
                    KM_crash_across_depth = squeeze(KM_crash(:, a, s)); 

                    % shields number exceed Komar threshold 
                    entrained_depth_index = find(shield_across_depth > KM_crash_across_depth, 1, 'first');

                    if isempty(entrained_depth_index)
                        d_crash{s, c}(u, a) = NaN;  % Does not reach entrainment
                    else
                        % Check if the wave breaks before entrainment happens
                        if H_shoal(entrained_depth_index) / L_shoal(entrained_depth_index) < max_steepness
                            d_crash{s, c}(u, a) = d(entrained_depth_index);  % Entrainment happens before breaking
                        else
                            d_crash{s, c}(u, a) = -999;  % Reaches entrainment after breaking occurs
                        end
                    end
                end
            end




    end
end

icecolor = {'#42F2F7','#46ACC2','#498C8A','#4B6858'}; % ice colors (blues)
organicicecolor = {'#ffa5ab','#da627d','#a53860','#450920'}; % organic-ice (reds)
organiccolor =  {'#f8ed62','#e9d700','#dab600','#a98600'}; % organics (yellows)

figure
for s = 2%1:numel(rho_s)

    if s == 1
        mycolor = organicicecolor;
    elseif s == 2
        mycolor = icecolor;
    elseif s == 3
        mycolor = organiccolor;
    end

    for c = 2%1:numel(lakes)

        depth_entrain_sand = d_crash{s,c}(:,1);
        depth_entrain_sand(depth_entrain_sand==-999) = NaN;
        depth_entrain_grav = d_crash{s,c}(:,end);
        depth_entrain_grav(depth_entrain_grav==-999) = NaN;
        plot(test_speeds,depth_entrain_sand,'--','Color','k','LineWidth',3,'DisplayName',[lakes{c} ' sand (rho = ' num2str(rho_s(s)) ')'])
        hold on
        plot(test_speeds,depth_entrain_grav,'-','Color','k','LineWidth',3,'DisplayName',[lakes{c} ' gravel(rho = ' num2str(rho_s(s)) ')'])
        sand1 = d_crash{1,c}(:,1)';
        sand1(sand1==-999) = NaN; 
        sand2 = d_crash{3,c}(:,1)';
        sand2(sand2==-999) = NaN;
        grav1 = d_crash{1,c}(:,end)';
        grav1(grav1==-999) = NaN;
        grav2 = d_crash{3,c}(:,end)';
        grav2(grav2==-999) = NaN;
        fill([test_speeds, fliplr(test_speeds)], [sand2, fliplr(sand1)],'k','HandleVisibility','off');
        fill([test_speeds, fliplr(test_speeds)], [grav2, fliplr(grav1)],'k','HandleVisibility','off');

    end

    legend('show','Location','best','Interpreter', 'latex')
    grid on;
    xlabel('wind speed [m/s]','FontSize', label_font_size,'Interpreter', 'latex')
    ylabel('entrainment depth [m]','FontSize', label_font_size,'Interpreter', 'latex')

end

for i = 1:numel(frac_time)
    xline(edges(i),'-k',{num2str((frac_time(i)))},'HandleVisibility','off','Interpreter', 'latex','FontSize',label_font_size)
end

avg_ol = [0.39 1.42 4.05];
avg_ol_label = {'summer daily avg','summer storm avg','summer max storm'};
for i = 1:numel(avg_ol)
    xline(avg_ol(i),'-r',avg_ol_label{i},'HandleVisibility','off','Interpreter', 'latex','FontSize',label_font_size,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
end

ws_sand = calc_settling_velocity(rho_s(2),d50(1),Planet);
ws_gravel = calc_settling_velocity(rho_s(2),d50(end),Planet);



figure;
plot(test_speeds,(depth_entrain_sand./(ws_sand))./(60*60),'--','Color','k','LineWidth',3)
hold on;
plot(test_speeds,(depth_entrain_grav./(ws_gravel))./(60*60),'-','Color','k','LineWidth',3)
ylabel('maximum settling time (hr)','FontSize', label_font_size,'Interpreter', 'latex')
xlabel('wind speed [m/s]','FontSize', label_font_size,'Interpreter', 'latex')
legend('fines','gravel')

for i = 1:numel(frac_time)
    xline(edges(i),'-k',{num2str((frac_time(i)))},'HandleVisibility','off','Interpreter', 'latex','FontSize',label_font_size)
end


for i = 1:numel(avg_ol)
    xline(avg_ol(i),'-r',avg_ol_label{i},'HandleVisibility','off','Interpreter', 'latex','FontSize',label_font_size,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
end

figure;
for i = 1:numel(lakes)
    plot(test_speeds,H0(i,:),'-','LineWidth',5,'Color',lakecolors(i,:),'DisplayName',lakes{i})
    hold on;
end
xlabel('wind speed [m/s]','FontSize', label_font_size)
ylabel('wave height [m]','FontSize', label_font_size)
grid on
legend('show')
for i = 1:numel(frac_time)
    xline(edges(i),'-k',{num2str((frac_time(i)))},'HandleVisibility','off','Interpreter', 'latex','FontSize',label_font_size)
end

for i = 1:numel(avg_ol)
    xline(avg_ol(i),'-r',avg_ol_label{i},'HandleVisibility','off','Interpreter', 'latex','FontSize',label_font_size,'LabelVerticalAlignment','middle','LabelHorizontalAlignment','center')
end



% WOLMAN-MILLER PLOT

u = round(mag_wind,1);
psi = round(angle_wind,1);

figHandle = figure;
yyaxis left
h = histogram(u,'BinEdges',0:0.1:max(u), 'Normalization','probability','FaceColor','none','LineStyle','--','LineWidth',1);
axisHandle = figHandle.Children;
histHandle = axisHandle.Children;
histHandle.BinEdges = histHandle.BinEdges + histHandle.BinWidth/2;
aa = histHandle.BinEdges;
percent_time = h.Values;
grid on
xlabel('wind speed [m/s]')
ylabel('PDF')

yyaxis right

activity = depth_entrain_sand; % 90 is max de
plot(test_speeds,activity,'LineWidth',5,'Color',lakecolors(2,:))
ylabel('Entrainment Depth [m]')
%ylim([0 1])

trimmed_fraction = percent_time;
trimmed_fraction(1) = [];
trimmed_fraction(end) = [];

% Wolman_Miller = activity.*trimmed_fraction';
% Wolman_Miller = Wolman_Miller / sum(Wolman_Miller,"omitnan");
nBins = length(Wolman_Miller);
x = (0:nBins-0.9).*0.1 + 0.1;
yyaxis left
hold on;
b1 = bar(x, Wolman_Miller,1,'FaceColor','r','FaceAlpha',0.3,'LineWidth',3);

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

function settling_velocity = calc_settling_velocity(rho_s,D50,Planet)


    % non-dimensional diameter
    Re_particle = (sqrt(((rho_s/Planet.rho_liquid)-1)*Planet.gravity*D50)*D50)/Planet.nu_liquid;
    D_star = Re_particle.^2;
    % non-dimensional settling velocity
    log_nondim_w_s = -3.76715 + 1.92944.*(log10(D_star)) - 0.09815.*(log10(D_star)).^2 - 0.00575.*(log10(D_star)).^3 + 0.00056.*(log10(D_star)).^4;
    nondim_w_s = 10.^(log_nondim_w_s);
    % dimensional settling velocity
    settling_velocity = ((((rho_s/Planet.rho_liquid)-1)*Planet.gravity*Planet.nu_liquid)*nondim_w_s)^(1/3);

end


% % 1. Define histogram data
% data = randn(1000,1);                         % some example data
% [hist_counts, edges] = histcounts(data, 50); % histogram with 50 bins
% 
% % 2. Compute bin centers from edges
% bin_centers = edges(1:end-1) + diff(edges)/2;
% 
% % 3. Define the curve y at bin centers
% % (e.g., a Gaussian curve or any function)
% y = exp(-bin_centers.^2); % example curve
% 
% % 4. Multiply the curve by histogram counts
% new_hist = hist_counts .* y;
% 
% % 5. Plot original histogram and the result
% figure;
% bar(bin_centers, hist_counts, 'FaceAlpha', 0.5); hold on;
% plot(bin_centers, y*max(hist_counts), 'r', 'LineWidth', 2); % scaled for visibility
% bar(bin_centers, new_hist, 'FaceAlpha', 0.5); % new histogram
% legend('Original Histogram', 'Curve (scaled)', 'New Histogram');
% xlabel('Bin Centers'); ylabel('Counts');
