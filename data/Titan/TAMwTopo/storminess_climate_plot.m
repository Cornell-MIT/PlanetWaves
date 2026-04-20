clc
clear
close all


all_lake = {'Ontario Lacus','Ligeia Mare'};   % lakes of interest 

%% --- Time Conversions --- %%
sec_per_Earth_day          = 86400;                                        % number of seconds in a single Earth day
Earth_days_per_year        = 365;                                          % number of Earth days in a single Earth year
TitanDay_EarthDays         = 15.945;                                       % number of Earth days in a single Titan day
TitanYear_EarthYears       = 29.5;                                         % number of Earth years in a single Titan year
TitanYear_EarthDays        = TitanYear_EarthYears * Earth_days_per_year;   % conversion bw Titan year and Earth year
Titan_days_per_Titan_year  = TitanYear_EarthDays / TitanDay_EarthDays;     % number of Titan days in a single Titan year
Hours_per_Titan_day        = TitanDay_EarthDays * 24;                      % number of hours in a Titan day

%% --- Modeling Sampling Info from Juan --- %%
years = 10;                                                                % data included 10 Titan years
sampling_freq = 1/6;                                                       % 1 sample per 6 hours
points_per_day = round(Hours_per_Titan_day * sampling_freq);
points_per_year = round(Titan_days_per_Titan_year * points_per_day);
total_points = years * points_per_year;
points_per_season = floor(points_per_year / 4);

get_original = 1;                                                          % 1 = pull original NC dataset which takes forever or 0 = load saved .mat

% pretty colors 
lake_colors = {'#DA344D','#3F88C5'};                                       % {Ontario Lacus, Ligeia Mare}
season_colors = {'#FF4D6A', '#4DA6FF', '#33FF99', '#FFE34D'};              % {autumn, winter, spring, summer}

fig_fft = figure('Name','FFT of Wind Magnitude'); hold on; grid on; box on;
xlabel('Period (Titan Days)'); ylabel('|FFT|');
set(gca,'XScale','log','YScale','log','FontSize',20);
%% --- Loop over lakes --- %%

for p = 1:numel(all_lake)
    lake = all_lake{p};
    
    fig_mag = figure('Name',['Wind Mag Time Series at ',lake]); ax_mag = axes; hold on; grid on; box on;
    xlabel('Time [Titan Days]'); ylabel('Wind speed');

    fig_dir = figure('Name',['Wind Dir Time Series at ',lake]); ax_dir = axes; hold on; grid on; box on;
    xlabel('Time [Titan Days]'); ylabel('Wind direction [rad]');

    fig_pol = figure('Name',['Polar Plot of Wind Time Series at ',lake]); ax_pol = polaraxes; hold on;

    %% --- Load Wind Data --- %%
    if get_original
        warning('This will take forever to load, but here you go...');
        %fn = "C:\Users\Owner\OneDrive\Documents\00_Main\Work\MIT\External_Data\TAM.10mwinds.L24.table950.k1e-4.y181-190.6-hourly.nc";
        fn = fullfile(pwd, '..', '..', '..', '..', '..', 'MIT', 'External_Data', 'TAM.10mwinds.L24.table950.k1e-4.y181-190.6-hourly.nc');
        fileinfo = ncinfo(fn);
        lon = ncread(fn, 'lon');
        lat = ncread(fn, 'lat');
        uc = ncread(fn, 'ucomp',[1 1 1],[length(lon) length(lat) stopi]);
        vc = ncread(fn, 'vcomp',[1 1 1],[length(lon) length(lat) stopi]);

        if strcmp(lake,'Ontario Lacus')
            loni = 34;
            lati = 3; 
        elseif strcmp(lake,'Ligeia Mare')
            loni = 45; 
            lati = 31; 
        else
            error('need to specify a lake to load with coordinates')
        end
        u_ol = squeeze(uc(loni,lati,:));
        v_ol = squeeze(vc(loni,lati,:));
        mag_wind = sqrt(u_ol.^2 + v_ol.^2);                                % m/s
        angle_wind = deg2rad(mod(rad2deg(atan2(v_ol,u_ol)) + 360,360));    % radians

    else
        if strcmp(lake,'Ontario Lacus')
            load('OL_winds.mat','mag_wind','angle_wind')
        elseif strcmp(lake,'Ligeia Mare')
            load('LM_winds.mat','mag_wind','angle_wind')
        else
            error('Unknown lake');
        end
    end
    
    %% --- Seasonal & Time Indices --- %%
    L = length(mag_wind); 
    for yr = 1:years
        year_start = (yr-1)*points_per_year + 1;
        year_end   = min(year_start + points_per_year - 1, L);
        if year_start > L, break; end
        year_idx = year_start:year_end;
        mag_year = mag_wind(year_idx);
        t_start = (yr-1)*points_per_year;
        t = t_start + (1:length(year_idx));

        points_in_year = length(mag_year);
        points_season = floor(points_in_year/4);
        for s = 1:4
            s_start = (s-1)*points_season + 1;
            s_end   = min(s*points_season,length(mag_year));
            idx = s_start:s_end;
            
            plot(ax_mag,t(idx)/points_per_day, mag_year(idx), '.','Color',season_colors{s});
            plot(ax_dir,t(idx)/points_per_day, deg2rad(angle_wind(idx)), '.','Color',season_colors{s});
            polarplot(ax_pol, deg2rad(angle_wind(idx)), mag_year(idx), '.', 'Color', season_colors{s});
        end
        
        % Draw vertical line at year end
        if length(mag_year) == points_per_year
            xline(ax_mag,t(end)/points_per_day);
            xline(ax_dir,t(end)/points_per_day);
        end
    end
    
    %% --- Wind Rose --- %%
    figure('Name',['Wind Rose - ',lake]);
    h = wind_rose(wrapTo360(angle_wind + 180), mag_wind);  
   
    % fix up plot to look better
    for k = 1:length(h)
        if isprop(h(k), 'LineWidth')
            h(k).LineWidth = 2;         
        end
        if isprop(h(k), 'Color')
            h(k).Color = 'k';  
        end
    end
    ax = gca;
    ax.FontSize = 50;
    ax.FontWeight = 'bold';
    ax.GridLineStyle = '--';
    ax.GridAlpha = 0.7;
    
    %% --- Compute FFT --- %%
    dt = 6*3600; % 6 hours sampling interval in seconds
    Fs = 1/dt;
    Y = fft(mag_wind - mean(mag_wind));
    P2 = abs(Y/L);
    P1 = P2(1:floor(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = Fs*(0:floor(L/2))/L; % frequency [Hz]
    period_seconds = 1 ./ f;
    period_Titan_days = period_seconds / (TitanDay_EarthDays*sec_per_Earth_day);
    
    
    figure(fig_fft); 
    loglog(period_Titan_days, P1, 'Color', lake_colors{p}, 'LineWidth',2, 'DisplayName', lake);
end

figure(fig_fft);
xline(0.5,':k','1/2 Titan Day','LineWidth',1,'HandleVisibility','off','LabelVerticalAlignment','middle','FontSize',18,'LabelHorizontalAlignment','center');
OneTitanYr =  TitanYear_EarthDays / TitanDay_EarthDays;
xline(OneTitanYr,':k','1 Titan Year','LineWidth',1,'HandleVisibility','off','LabelVerticalAlignment','middle','FontSize',18,'LabelHorizontalAlignment','center');
legend('show','Location','southeast','FontSize',18);
xlim([1e-2 1e3]); % adjust period range
grid on;



