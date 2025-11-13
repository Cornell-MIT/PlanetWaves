clc
clear
close all

addpath(fullfile('..','past_runs'))
addpath(fullfile('..','..','planetwaves'))
addpath(fullfile('..','..','planetwaves','post_analysis'))

addpath(fullfile('..','..','data/Titan/TAMwTopo/'))
addpath(fullfile('..','..','data/Titan/TitanLakes/Bathymetries/bathtub_bathy'))


% Results from Titan_Waves03_43
load('../past_runs/Waves03_43.mat','test_speeds','time_to_run','wind_direction','zDep','buoy_loc','lakes','Planet','Model','H0','L0','T','C0','Cg0')
lakecolors = {'#F2D1C9','#E086D3','#8332AC','#462749'};

figure;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT ONE: PLOT OF ENTRAINMENT DEPTH VS WIND SPEED (w LAKE COMPOSITION)

icecolor = {'#42F2F7','#46ACC2','#498C8A','#4B6858'}; % ice colors (blues)
organicicecolor = {'#ffa5ab','#da627d','#a53860','#450920'}; % organic-ice (reds)
organiccolor =  {'#f8ed62','#e9d700','#dab600','#a98600'}; % organics (yellows)

figure
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
        plot(test_speeds,depth_entrain_sand,'-^','Color',mycolor{c},'LineWidth',3,'MarkerFaceColor',mycolor{c},'DisplayName',[lakes{c} ' sand (rho = ' num2str(rho_s(s)) ')'])
        hold on
        plot(test_speeds,depth_entrain_grav,'-s','Color',mycolor{c},'LineWidth',3,'MarkerFaceColor',mycolor{c},'DisplayName',[lakes{c} ' gravel(rho = ' num2str(rho_s(s)) ')'])
    end

    legend('show','Location','best')
    grid on;
    xlabel('wind speed [m/s]')
    ylabel('entrainment depth [m]')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT TWO: ENTRAINMENT DEPTH VS GRAIN SIZE (w LAKE COMPOSITION)


figure;
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
%% PLOT THREE: PLOT OF WOLMAN-MILLER GEOMORPHOLOGICAL EFFECTIVENESS 

% JUAN'S TAM WIND MODEL WINDS AT ONTARIO LACUS
% load('TAM_OL_winds.mat')
% str = '#F0C808';
% wind_color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
% 
% figHandle = figure;
% yyaxis left
% h = histogram(mag_wind,'BinEdges',0:0.1:max(mag_wind), 'Normalization','probability','FaceColor',wind_color);
% axisHandle = figHandle.Children;
% histHandle = axisHandle.Children;
% histHandle.BinEdges = histHandle.BinEdges + histHandle.BinWidth/2;
% aa = histHandle.BinEdges;
% percent_time = h.Values;
% grid on
% xlabel('wind speed [m/s]')
% ylabel('pdf time')
% 
% tt = [0 0 test_speeds];
% Qs = NaN(numel(lakes),numel(tt));
% max_angle = (cos(deg2rad(45))^(6/5))*sin(deg2rad(45));
% for c = [4 ol 1]
% 
%     yyaxis left
%     Hw = [0 0 H0(c,:)];
%     Tw = [0 0 T(c,:)];
%     Qs(c,:) = max_angle.*(Hw.^(5/6)).*(Tw.^(1/5));
%     Qs(c,:) = (Qs(c,:)./max(Qs(c,:)));
% 
%     Wolman_Miller(c,:) = Qs(c,:).*(percent_time*100);
% 
%     nBins = length(Wolman_Miller(c,:));
%     x = (0:nBins-0.9).*0.1 + 0.1;
%     yyaxis right
%     hold on;
%     if c == ol
%         b1 = bar(x, Wolman_Miller(c,:),1,'FaceColor',lakecolors{c},'FaceAlpha',0.3);
%     else
%         b1 = bar(x, Wolman_Miller(c,:),1,'FaceColor',lakecolors{c},'FaceAlpha',0,'LineStyle','--');
%     end
% 
% 
%     start_plot = find(Qs(c,:)==0,1,'last');
% 
%     plot(tt(start_plot:end),Qs(c,start_plot:end),':','LineWidth',2,'Color',lakecolors{c})
% 
% end
% ylabel('Qs')

% JUAN'S TAM WIND MODEL WINDS AT LIGEIA MARE
load('TAM_LM_winds.mat')
str = '#F0C808';
wind_color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

figHandle = figure;
yyaxis left
h = histogram(mag_wind,'BinEdges',0:0.1:4.5, 'Normalization','probability','FaceColor',wind_color);
axisHandle = figHandle.Children;
histHandle = axisHandle.Children;
histHandle.BinEdges = histHandle.BinEdges + histHandle.BinWidth/2;
aa = histHandle.BinEdges;
percent_time = h.Values;
grid on
xlabel('wind speed [m/s]')
ylabel('pdf time')

tt = [0 0 test_speeds];
Qs = NaN(numel(lakes),numel(tt));
max_angle = (cos(deg2rad(45))^(6/5))*sin(deg2rad(45));
for c = [4 lm 1]

    yyaxis left
    Hw = [0 0 H0(c,:)];
    Tw = [0 0 T(c,:)];
    Qs(c,:) = max_angle.*(Hw.^(5/6)).*(Tw.^(1/5));
    Qs(c,:) = (Qs(c,:)./max(Qs(c,:)));
    
    Wolman_Miller(c,:) = Qs(c,:).*(percent_time*100);

    nBins = length(Wolman_Miller(c,:));
    x = (0:nBins-0.9).*0.1 + 0.1;
    yyaxis right
    hold on;
    if c == lm
        b1 = bar(x, Wolman_Miller(c,:),1,'FaceColor',lakecolors{c},'FaceAlpha',0.3);
    else
        b1 = bar(x, Wolman_Miller(c,:),1,'FaceColor',lakecolors{c},'FaceAlpha',0);
    end

   
    start_plot = find(Qs(c,:)==0,1,'last');
    
    plot(tt(start_plot:end),Qs(c,start_plot:end),':','LineWidth',2,'Color',lakecolors{c})

end
ylabel('Qs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT FOUR: CONTOURS OF ENTRAINMENT IN ONTARIO LACUS
load('ol_bathtub_0.002000_slope.mat')

zDep(isnan(zDep))=-1;
zDep_entrain = zDep;
zDep_entrain(zDep_entrain>d_crash{ice,ol}(18,1)) = -1;
figure;
M2 = contourf(zDep_entrain,[0 d_crash{ice,ol}(18,1)],'-w','LineWidth',3);
colormap([0.5 0.5 0.5])
hold on;
M1 = contour(zDep,[0:10:max(max(zDep))],'-k','LineWidth',3,'ShowText','on');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
%exportgraphics(gcf, 'ol_smallsand_medwind.pdf', 'ContentType', 'vector');

load('lm_bathtub_0.002000_slope.mat')
zDep(isnan(zDep))=-1;
zDep_entrain = zDep;
zDep_entrain(zDep_entrain>d_crash{ice,lm}(18,1)) = -1;
figure;
M2 = contourf(zDep_entrain,[0 d_crash{ice,lm}(18,1)],'-w','LineWidth',3);
colormap([0.5 0.5 0.5])
hold on;
M1 = contour(zDep,[0:20:max(max(zDep))],'-k','LineWidth',3,'ShowText','on');

set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
%exportgraphics(gcf, 'lm_smallsand_maxwind.pdf', 'ContentType', 'vector');
