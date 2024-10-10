clc
clear
close all

addpath(fullfile('..','data/Titan/TAMwTopo/'))
addpath(fullfile('..','planetwaves'))
addpath(fullfile('..','planetwaves','post_analysis'))
addpath(fullfile('.','past_runs'))

% Results from Titan_Ontario_Lacus_Entrainment
load('Waves03_43.mat','H0','test_speeds','lakes','lake_slope','L0','T','Planet','C0','Cg0','rho_s','d50','ice','lakecolors','time_to_run','wind_direction','zDep','buoy_loc','max_steepness','min_depth')

alpha_0 = 0;
d = 100:-lake_slope:min_depth;
d50 = d50(1:100:end);
d50 = [d50(1) d50(end)];


for c = 1:numel(lakes)

    [Planet, Model, Wind, Uniflow, Etc] = initalize_model(lakes{c}, time_to_run, wind_direction, zDep, buoy_loc);

    % smooth artifacts in wave model (introduces max of ~ 0.1 error)
    % figure
    % plot(test_speeds,H0(1,:),'-','DisplayName','original H')
    % hold on
    H0(c,:) = smooth_waves(H0(c,:));
    % plot(test_speeds,H0(1,:),'--','DisplayName','smoothed H')
    % figure;
    % plot(test_speeds,T(1,:),'-','DisplayName','original T')
    % hold on
    T(c,:) = smooth_waves(T(c,:));
    % plot(test_speeds,T(1,:),'--','DisplayName','smoothed T')
    % figure;
    % plot(test_speeds,L0(1,:),'-','DisplayName','original L0')
    L0(c,:) = smooth_waves(L0(c,:));
    % hold on
    % plot(test_speeds,L0(1,:),'--','DisplayName','smoothed L0')


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
            um_shoal = (H_shoal / 2) .* ((Planet.gravity * T(c,u)) ./ L_shoal) .* (1 ./ cosh(2*pi*d_L));
            break_frac = H_shoal ./ L_shoal;
            break_frac(break_frac >= max_steepness) = NaN;

             % meshgrid to vectorizing for computational speed
            [S, A] = meshgrid(rho_s, d50);  % [S,A] : [numel(d50), numel(rho_s)]

            shields = zeros(length(um_shoal), size(S, 1), size(S, 2));  % size : [numel(d), numel(d50), numel(rho_s)]
            KM_crash = zeros(length(um_shoal), size(A, 1), size(A, 2));  % size : [numel(d), numel(d50), numel(rho_s)]

            for z = 1:length(um_shoal)
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
        plot(test_speeds,depth_entrain_sand,'^','Color',mycolor{c},'LineWidth',3,'MarkerFaceColor',mycolor{c},'DisplayName',[lakes{c} ' sand (rho = ' num2str(rho_s(s)) ')'])
        hold on
        plot(test_speeds,depth_entrain_grav,'s','Color',mycolor{c},'LineWidth',3,'MarkerFaceColor',mycolor{c},'DisplayName',[lakes{c} ' gravel(rho = ' num2str(rho_s(s)) ')'])
    end

    legend('show','Location','best')
    grid on;
    xlabel('wind speed [m/s]')
    ylabel('entrainment depth [m]')

end


% JUAN'S TAM WIND MODEL WINDS
load('TAM_winds.mat')
str = '#F0C808';
wind_color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

figHandle = figure;
yyaxis left
h = histogram(mag_wind,'BinEdges',0:0.1:max(mag_wind), 'Normalization','probability','FaceColor',wind_color);
axisHandle = figHandle.Children;
histHandle = axisHandle.Children;
histHandle.BinEdges = histHandle.BinEdges + histHandle.BinWidth/2;
aa = histHandle.BinEdges;
percent_time = h.Values;
grid on
xlabel('wind speed [m/s]')
ylabel('pdf time')

tt = [0 0 test_speeds];
Qc = NaN(numel(lakes),numel(tt));
for c = [4 3 2 1]

    yyaxis left
    Hw = [0 0 H0(c,:)];
    Tw = [0 0 T(c,:)];
    Qs(c,:) = 0.5.*(Hw.^(5/6)).*(Tw.^(1/5));
    Qs(c,:) = (Qs(c,:)./max(Qs(c,:)));
    
    Wolman_Miller(c,:) = Qs(c,:).*(percent_time*100);

    nBins = length(Wolman_Miller(c,:));
    x = (0:nBins-0.9).*0.1 + 0.1;
    yyaxis right
    hold on;
    b1 = bar(x, Wolman_Miller(c,:),1,'FaceColor',lakecolors{c});
    b1.FaceAlpha = 0.3;
    % plot(x,Wolman_Miller(c,:),'--','Color',lakecolors{c},'LineWidth',2)
    
    
   
    start_plot = find(Qs(c,:)==0,1,'last');
    
    plot(tt(start_plot:end),Qs(c,start_plot:end),':','LineWidth',2,'Color',lakecolors{c})

end
ylabel('Qs')

