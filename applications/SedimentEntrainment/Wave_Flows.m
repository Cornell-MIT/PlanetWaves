clc
clear
close all


disp('====== RUNNING WAVE FLOW ENTRAINMENT ==================')
% calculates sediment entrainment under shoaling waves at Titan

run_waves = 0; % if true, will run PlanetWaves (Titan_DeepwaterWaves.m), otherwise will load wave results from past run

% wave model
addpath(genpath(fullfile('..','..','planetwaves')))
% GCM surface winds
addpath(fullfile('..','..','data/Titan/TAMwTopo/'))
% simple bathymetry for Ontario Lacus using constant slope
addpath(fullfile('..','..','data/Titan/TitanLakes/Bathymetries/bathtub_bathy'))

if run_waves
    Titan_DeepwaterWaves;
else
    pathname = [pwd,'\past_runs\Titan_DeepWaterWaves.mat'];
    last_created = dir(pathname).date;
    fprintf('Loading data from file run: %s\n',last_created);
    % Results from last run of Titan_DeepwaterWaves
    load('./past_runs/Titan_DeepwaterWaves.mat','test_speeds','time_to_run','wind_direction','zDep','buoy_loc','lakes','Planet','Model','H0','L0','T','C0','Cg0')
end

lakecolors = {'#F2D1C9','#E086D3','#8332AC','#462749'};

figure('Name','u v H inputs');
for i = 1:numel(lakecolors)
    plot(test_speeds,H0(i,:),'-','LineWidth',5,'Color',lakecolors{i},'DisplayName',lakes{i})
    hold on;
end
xlabel('wind speed [m/s]')
ylabel('wave height [m]')
grid on
legend('show')


lake_slope = 2e-3;
alpha_0 = 0;
min_depth = 10^-4;
d50 = 6.35e-5:1e-5:0.1;  % [fines gravel]

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
    H0(c,:) = smooth_waves(H0(c,:));
    T(c,:) = smooth_waves(T(c,:));
    L0(c,:) = smooth_waves(L0(c,:));

    for u = 1:numel(test_speeds)

            % shoaling waves
            L_shoal = L0(c,u) * sqrt(tanh(4*(pi^2)*d/(T(c,u)^2*Planet.gravity))); % Eckhart aproximation to not solve recursively 
            d_L = d ./ L_shoal;
      
            alpha = asin(tanh((2*pi.*d_L).*sin(alpha_0)));
            C_shoal = C0(c,u) * tanh(2*pi*d_L);
            n = 0.5 * (1 + ((4*pi*d_L)./sinh(4*pi*d_L)));
            Cg_shoal = n .* C_shoal;
            KR = sqrt(cos(alpha_0)./cos(alpha));
            KS = sqrt(Cg0(c,u)./Cg_shoal);
            H_shoal = KR .* KS .* H0(c,u);
            % orbital size, velocity, breaking wave condition
            d0_shoal = (H_shoal / 2) .* (1 ./ sinh(2*pi*d_L)); % <-- from Komar/Miller, for z = -d
            if test_speeds(u) == 0.7 || test_speeds(u) == 0.8
                max_d0_shoal = [max_d0_shoal max(d0_shoal,[],'omitnan')];
            end
            um_shoal = (H_shoal ) .* ((Planet.gravity * T(c,u)) ./ L_shoal) .* (1 ./ cosh(2*pi*d_L)); % <-- from Komar/Miller, deep-water
            break_frac = H_shoal ./ L_shoal;
            break_frac(break_frac >= max_steepness) = NaN;

             % meshgrid to vectorizing for computational speed
            [S, A] = meshgrid(rho_s, d50);  % [S,A] : [numel(d50), numel(rho_s)]

            shields = zeros(length(d), size(S, 1), size(S, 2));  % size : [numel(d), numel(d50), numel(rho_s)]
            KM_crash = zeros(length(d), size(A, 1), size(A, 2));  % size : [numel(d), numel(d50), numel(rho_s)]
            KM_turb_crash =  zeros(length(d), size(A, 1), size(A, 2));
            for z = 1:length(d)
                shields(z, :, :) = (Planet.rho_liquid * (um_shoal(z)^2)) ./ ((S - Planet.rho_liquid) * Planet.gravity .* A);
                KM_crash(z, :, :) =   0.341300 .* sqrt(d0_shoal(z) ./ A);    % fit for laminar data in Figure 4, Komar & Miller 1973 (see boundary_layer_threshold.m in estimating boundary layer folder)
                KM_turb_crash(z, :, :) = 0.182849 .* sqrt(d0_shoal(z) ./ A); % fit for laminar data in Figure 4, Komar & Miller 1973 (see boundary_layer_threshold.m in estimating boundary layer folder)
            end

            for s = 1:numel(rho_s)

                for a = 1:numel(d50)

                    % across all depths                
                    shield_across_depth = squeeze(shields(:, a, s));  
                    KM_crash_across_depth = squeeze(KM_crash(:, a, s)); 
                    KM_turb_crash_across_depth = squeeze(KM_turb_crash(:,a,s));

                    
                    % shields number exceed Komar threshold 
                    entrained_depth_index_laminar = find(shield_across_depth > KM_crash_across_depth, 1, 'first');
                    entrained_depth_index_turb = find(shield_across_depth > KM_turb_crash_across_depth, 1, 'first');

                    % figure('Name','Shields Curves');
                    % plot(d,shield_across_depth,'-k')
                    % hold on
                    % plot(d,KM_crash_across_depth,'-r')
                    % plot(d,KM_turb_crash_across_depth,'-b')
                    % plot(d(entrained_depth_index_laminar),shield_across_depth(entrained_depth_index_laminar),'or')
                    % plot(d(entrained_depth_index_turb),shield_across_depth(entrained_depth_index_turb),'ob')
                    % set(gca,'YScale','log')
                    % ylabel('$\rho * u_m^2 / ((\rho_s - \rho)gD_{50}$','Interpreter','latex')
                    % xlabel('depth')
                    
                    if ~isempty(entrained_depth_index_turb) && ~isempty(entrained_depth_index_laminar)
                        Re_particle_opt1 = (um_shoal(entrained_depth_index_laminar)*d50(a))/Planet.nu_liquid;
                        Re_particle_opt2 = (um_shoal(entrained_depth_index_turb)*d50(a))/Planet.nu_liquid;
    
                        % if particle Re says boundary layer is turbulent, use turbulent boundary prediction curve
                        % else, use laminar boundary prediction curve
                        if Re_particle_opt2 > 100 && Re_particle_opt1 > 100
                            entrained_depth_index = entrained_depth_index_turb;
                        elseif Re_particle_opt1 > 100 
                            entrained_depth_index = entrained_depth_index_turb;
                        else
                            entrained_depth_index = entrained_depth_index_laminar;
                        end
                    elseif ~isempty(entrained_depth_index_turb)
                        
                        entrained_depth_index = entrained_depth_index_turb;

                    elseif ~isempty(entrained_depth_index_laminar)

                        entrained_depth_index = entrained_depth_index_laminar;    

                    else
                        entrained_depth_index = [];
                    end
  
                    if isempty(entrained_depth_index)
                        d_crash{s, c}(u, a) = NaN;  % Does not reach entrainment
                        boundary_layer_RE{s, c}(u, a) = NaN;
                    else
                        % Check if the wave breaks before entrainment happens
                        if H_shoal(entrained_depth_index) / L_shoal(entrained_depth_index) < max_steepness
                            d_crash{s, c}(u, a) = d(entrained_depth_index);  % Entrainment happens before breaking
                            boundary_layer_RE{s, c}(u, a) = (um_shoal(entrained_depth_index)*d50(a))/Planet.nu_liquid;
                        else
                            d_crash{s, c}(u, a) = -999;  % Reaches entrainment after breaking occurs
                            boundary_layer_RE{s, c}(u, a) = -999;
                        end
                    end
                end
            end




    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT ONE: PLOT OF ENTRAINMENT DEPTH VS WIND SPEED (w LAKE COMPOSITION)

icecolor = {'#42F2F7','#46ACC2','#498C8A','#4B6858'}; % ice colors (blues)
organicicecolor = {'#ffa5ab','#da627d','#a53860','#450920'}; % organic-ice (reds)
organiccolor =  {'#f8ed62','#e9d700','#dab600','#a98600'}; % organics (yellows)

figure('Name','u vs entrainment depth')
hold on
xline(0.39,'-','summer daily average','LabelHorizontalAlignment','left','LineWidth',4,'Color','k','HandleVisibility','off','LabelVerticalAlignment','middle','FontSize',20)
xline(1.42,'-','summer storm average','LabelHorizontalAlignment','left','LineWidth',4,'Color','k','HandleVisibility','off','LabelVerticalAlignment','middle','FontSize',20)
xline(4.05,'-','summer maximum','LabelHorizontalAlignment','left','LineWidth',4,'Color','k','HandleVisibility','off','LabelVerticalAlignment','middle','FontSize',20)
for s = 1:numel(rho_s)

    if s == 1
        mycolor = organicicecolor;
    elseif s == 2
        mycolor = icecolor;
    elseif s == 3
        mycolor = organiccolor;
    end

    for c = 1:numel(lakes)

        depth_entrain_sand = d_crash{s,c}(:,1);
        depth_entrain_sand(depth_entrain_sand==-999) = NaN;
        depth_entrain_grav = d_crash{s,c}(:,end);
        depth_entrain_grav(depth_entrain_grav==-999) = NaN;

        compare_sand{s}(c,:) = depth_entrain_sand;
        compare_grav{s}(c,:) = depth_entrain_grav;
        
        plot(test_speeds,depth_entrain_sand,':','Color',mycolor{c},'LineWidth',4,'MarkerFaceColor',mycolor{c},'DisplayName',[lakes{c} ' sand (rho = ' num2str(rho_s(s)) ')'])
        plot(test_speeds,depth_entrain_grav,'-','Color',mycolor{c},'LineWidth',4,'MarkerFaceColor',mycolor{c},'DisplayName',[lakes{c} ' gravel(rho = ' num2str(rho_s(s)) ')'])
    end

    legend('show','Location','best')
    grid on;
    xlabel('Wind speed [m/s]')
    ylabel('Maximum entrainment depth [m]')
    xlim([0 4.5])

end

box on;
ax = gca;
ax.FontSize = 16;

avg_sand_diff = mean(compare_sand{2}(2,:)-compare_sand{2}(3,:),'omitnan');
avg_grav_diff = mean(compare_grav{2}(2,:)-compare_grav{2}(3,:),'omitnan');
std_sand_diff = std(compare_sand{2}(2,:)-compare_sand{2}(3,:),'omitnan');
std_grav_diff = std(compare_grav{2}(2,:)-compare_grav{2}(3,:),'omitnan');

fprintf('Average difference:\n%f pm %f m (sand)\n%f pm %f (gravel)\n',avg_sand_diff,std_sand_diff,avg_grav_diff,std_grav_diff)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT TWO: ENTRAINMENT DEPTH VS GRAIN SIZE (w LAKE COMPOSITION)


figure('Name','d50 vs entrainment depth');
for c = 1:numel(lakes)
    entrainment_depth = d_crash{ice,c}(end,:);
    entrainment_depth(entrainment_depth==-999) = NaN;
    plot(d50,d_crash{ice,c}(end,:),'-','LineWidth',3,'Color',lakecolors{c},'DisplayName',lakes{c})
    hold on
end
legend('show')
xlabel('d50')
ylabel('entrainment depth')
grid on;
hold off
set(gca,'XScale','log')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOT 3: ENTRAINMENT DEPTH FOR SMALL ICE GRAINS IN ONTARIO LACUS
[~,summer_avg] = min(abs(test_speeds - 1.42)); % use summer wind average (~1.42 m/s)
ice_entrain_summer_average = d_crash{ice,2}(summer_avg,1);
map_entrainment_depth_in_Ontario_Lacus(ice_entrain_summer_average)

fprintf('Ice fine sand grains in Ontario Lacus composition %0.1f m\n',ice_entrain_summer_average)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT 4: Entrainment Frequency PLot

% climate model results for surface wind speeds
load('OL_winds','mag_wind','angle_wind');

% Isolate southern summer winds
points_per_day = 64;
titan_days_per_year = 674;         
points_per_year = titan_days_per_year * points_per_day;
points_per_season = floor(points_per_year / 4);
years = floor(length(mag_wind) / points_per_year);
autumn_idx = cell(years,1);
winter_idx = cell(years,1);
spring_idx = cell(years,1);
summer_idx = cell(years,1);
for y = 1:years
    year_start = (y-1)*points_per_year + 1;
    autumn_idx{y} = year_start : year_start + points_per_season - 1;
    winter_idx{y} = year_start + points_per_season : year_start + 2*points_per_season - 1;
    spring_idx{y} = year_start + 2*points_per_season : year_start + 3*points_per_season - 1;
    summer_idx{y} = year_start + 3*points_per_season : year_start + 4*points_per_season - 1;
end
autumn_all = [autumn_idx{:}];
winter_all = [winter_idx{:}];
spring_all = [spring_idx{:}];
summer_all = [summer_idx{:}];
wind_S_summer = mag_wind(winter_all); % northern winter = southern summer

% Calculate fraction of time for max entrainment
figure('Name','Rarity of observation');
subplot(2,2,[2 4])
hold on
grid on
ylabel('Minimum entrainment depth [m]', 'Interpreter', 'latex', 'FontSize',25)
xlabel('Fraction of time (\%)','Interpreter','latex','FontSize',25)
set(gca,'XScale','log')
xlim([10^-2 200])
ylim([0 11])
set(gca,'YDir','reverse')
yline(ice_entrain_summer_average,'-','summer storm average (u$_{10}$ = 1.42 m/s)','Color',[0.5 0.5 0.5],'LineWidth',2,'HandleVisibility','off','Interpreter','latex')
ax = gca;
ax.FontSize = 20; 

wind_label = {'All Seasons','Southern Summer'};
for wind_climate = 1:2 % first loop is all winds, second is just the subset for southern summer

    %  Entrainment Depth vs Fraction of time
    if wind_climate == 1
        wind = mag_wind;
        linecolor = 'k';
    elseif wind_climate == 2
        wind = wind_S_summer;
        linecolor = 'r';
    end

    % make wind speed bins 
    bin_size = 0.1; % m/s
    edges_hist = 0.0:bin_size:max(wind)+bin_size;                           % edge of bins for histcount 
    [counts, ~] = histcounts(wind, edges_hist);                             % count up number of occurences within each bin
    

    frac_time_bins = counts / sum(counts);                                  % fraction of time per bin
    bin_centers = edges_hist(1:end-1) + diff(edges_hist)/2;                 % assign wind speed to value at center of bin
    

    depth_entrain_grav(isnan(depth_entrain_grav)) = 0;
    depth_entrain_sand(isnan(depth_entrain_sand)) = 0;
    depth_entrain_sand_interp  = interp1(test_speeds, depth_entrain_sand, bin_centers, 'linear', NaN);
    depth_entrain_grav_interp = interp1(test_speeds, depth_entrain_grav, bin_centers, 'linear', NaN);
    depth_entrain_sand_interp(isnan(depth_entrain_sand_interp)) = 0;
    depth_entrain_grav_interp(isnan(depth_entrain_grav_interp)) = 0;

    %cumulative time at least depth Y 
    cum_time_sand = flip(cumsum(flip(frac_time_bins)));
    cum_time_grav = cum_time_sand;  % same time distribution
    % Add zero-depth no-move time to beginning
    xvals = [1; cum_time_sand'];   % start at 100% at zero depth
    y_sand = [0; depth_entrain_sand_interp'];
    y_grav = [0; depth_entrain_grav_interp'];
    

    plot(xvals.*100, y_sand, ':', 'LineWidth', 3, 'Color', linecolor, 'DisplayName',['Fine Sand (63.5 $\mu$m), ', wind_label{wind_climate}]);
    plot(xvals.*100, y_grav, '-', 'LineWidth', 3, 'Color', linecolor,'DisplayName',['Cobbles (0.1 m), ', wind_label{wind_climate}]);

    [y_sand_unique,i_unique,~] = unique(y_sand);
    x_vals_sand_unique = xvals(i_unique).*100;
    [y_grav_unique,i_unique,~] = unique(y_grav);
    x_vals_grav_unique = xvals(i_unique).*100;

    time_at_1m_sand(wind_climate) = interp1(y_sand_unique,x_vals_sand_unique,1);
    time_at_1m_grav(wind_climate) = interp1(y_grav_unique,x_vals_grav_unique,1);
    time_at_6p6m_sand(wind_climate) = interp1(y_sand_unique,x_vals_sand_unique,6.6);
    time_at_6p6m_grav(wind_climate) = interp1(y_grav_unique,x_vals_grav_unique,6.6);
    
    
    
end

box on;
legend('show','Location','best','Interpreter','latex')
hold off

disp('===========')
fprintf('Sand time at 1-m depth: %0.1f (all) -- %0.1f (southern summer) percent of time\n',time_at_1m_sand)
fprintf('Cobble time at 1-m depth: %0.1f (all) -- %0.1f (southern summer) percent of time\n',time_at_1m_grav)
fprintf('Sand time at 6.6-m depth: %0.1f (all) -- %0.1f (southern summer) percent of time\n',time_at_6p6m_sand)
fprintf('Cobble time at 6.6-m depth: %0.1f (all) -- %0.1f (southern summer) percent of time\n',time_at_6p6m_grav)
disp('===========')

% wind speed histogram
edges = 0:0.01:max(mag_wind);
centers = edges(1:end-1) + diff(edges)/2;
no_waves = mag_wind < 0.5;
no_waves_count = histcounts(mag_wind(no_waves), edges)/ numel(mag_wind);
waves_count = histcounts(mag_wind(~no_waves), edges)/ numel(mag_wind);
subplot(2,2,1)
stairs(centers, no_waves_count, 'k', 'LineWidth', 1.5)
hold on
bar(centers, waves_count, 1, 'FaceColor','k','EdgeColor','k')
xlabel('wind speed, u$_{10}$ [m/s]','Interpreter','latex','FontSize',25)
ylabel('Fraction of time','Interpreter','latex','FontSize',25)
xlim([0 5])
grid on;
box on;
ax = gca;
ax.FontSize = 20; 

% wind speed vs wave height
subplot(2,2,3)
for i = 1:numel(lakecolors)
    H = H0(i,:);
    H(H==0) = NaN;
    plot(test_speeds,H,'-','LineWidth',5,'Color',lakecolors{i},'DisplayName',lakes{i})
    hold on;
end
xline(0.5,':k','HandleVisibility','off')
xlabel('wind speed, u$_{10}$ [m/s]','Interpreter','latex','FontSize',25)
ylabel('wave height [m]','Interpreter','latex','FontSize',25)
grid on
legend('show')
xlim([0 5])
box on;
ax = gca;
ax.FontSize = 20; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
% % PLOT 5: Boundary Layer Reynold's Number (Fig 6 Komar & Miller 1973)
% 
% speed_color = winter(numel(test_speeds));
% figure('Name','D vs Re')
% hold on
% yline(100)
% for s = 2%1:numel(rho_s)
% 
%     if s == 1
%         mycolor = organicicecolor;
%     elseif s == 2
%         mycolor = icecolor;
%     elseif s == 3
%         mycolor = organiccolor;
%     end
% 
%     for c = 2%1:numel(lakes)
% 
%         for u = 1:numel(test_speeds)
%             RE_waves = boundary_layer_RE{s, c}(u, :);
%             plot(d50,RE_waves,'Color',speed_color(u,:));
%         end
% 
%     end
% 
% end
% xlabel('D50')
% ylabel('Re = u_mD/v')
% set(gca,'YScale','log')
% set(gca,'XScale','log')