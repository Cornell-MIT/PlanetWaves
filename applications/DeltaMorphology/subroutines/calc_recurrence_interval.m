function [reoccur_interval,x_sorted] = calc_recurrence_interval(x_time_series, x_target)
% Calculates the reoccurence interval for wind speeds with an even (6-hour) spacing interval
% in terms of Titan years. Reoccurence interval is the amount of time (in years) divided
% by the number of events of a given magnitude.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%   x_time_series : wind time series (assuming 6-hour evenly spaced intervals as seen in TAM [m/s]
%   x_target      : wind speed that I want to find the reocurrence interval for (e.g., 3.7) [m/s]
% OUTPUTS:
%   reoccur_interval : reoccurence interval [Titan years] 
%   x_sorted         : sorted list of annual maximums [m/s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example call: 
% [r,x_sorted] = calc_recurrence_interval(mag_wind,3.7) will report the reoccurence interval 
% for 3.7 in the time series mag_wind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    x_time_series = x_time_series(:);


    % time conversion to get in terms of Titan years
    Earth_days_per_year        = 365;                                          % number of Earth days in a single Earth year
    TitanDay_EarthDays         = 15.945;                                       % number of Earth days in a single Titan day
    TitanYear_EarthYears       = 29.5;                                         % number of Earth years in a single Titan year
    TitanYear_EarthDays        = TitanYear_EarthYears * Earth_days_per_year;   % conversion bw Titan year and Earth year
    Titan_days_per_Titan_year  = TitanYear_EarthDays / TitanDay_EarthDays;     % number of Titan days in a single Titan year
    Hours_per_Titan_day        = TitanDay_EarthDays * 24;                      % number of hours in a Titan day
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TAM SPECIFIC PARAMETERS (not an input for now)
    sampling_freq = 1/6;                                                       % 1 sample per 6 hours
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    points_per_day = round(Hours_per_Titan_day * sampling_freq);
    points_per_year = round(Titan_days_per_Titan_year * points_per_day);
    total_points = points_per_year;


    % how many years in data (trim remainders)
    n_total = length(x_time_series);
    n_years = floor(n_total / total_points);

    if n_years < 2
        error('Not enough data for recurrence analysis. Need more than 2 years.');
    end

    % remove any stragglers for incomplete years 
    x_time_series = x_time_series(1:n_years * points_per_year);
    
    number_above_target = sum(x_time_series>=x_target);
    fprintf('There are %i events with magnitude greater than %.2f over %i years\n',number_above_target,x_target,n_years);

    % matrix where columns are years and rows are the values
    % e.g.,
    %  yearly_matrix    = [ 10 5 9 7
    %                       1  7 4 1
    %                       2  1 0 1]
    % would describe four years where the third data point for each year
    % is 2 1 0 1
    yearly_matrix = reshape(x_time_series, points_per_year, n_years); 
    % Calculate maximum for each year
    % e.g.,
    % annual_max  = [10 7 9 7] from yearly_matrix
    annual_max = max(yearly_matrix, [], 1)';
    
    figure;
    plot(1:numel(annual_max),annual_max,'-ok','LineWidth', 1.5)
    xlabel('year [Titan years]')
    ylabel('annual max wind speed [m/s]')
    ylim([0 max(annual_max+1)])
    grid on;
    yline(x_target,'--r')

    [counts,bins] = hist(x_time_series(x_time_series>min(annual_max)));
    figure; 
    barh(bins,counts)
    xlabel('number exceeding minimum annual max')
    ylabel('wind speed')
    % e.g., 
    % x_sorted    = [10 9 7 7] 
    [x_sorted,rank] = sort(annual_max, 'descend'); % biggest has a rank of 1
    % Calculate the reoccurence interval for sorted and ranked annual maximum
    P = rank ./ (n_years + 1); % (probability) Weibull plotting position 

    figure;
    plot(P,x_sorted,'-k','LineWidth', 1.5)
    ylabel('annual max wind speed')
    xlabel('exceedance probability')

    reoccur_interval = 1./P;
    fprintf('------------------------------------------------------------------------------------------------\n')
    fprintf('Year  \t    Annual Max  \t   Rank  \t Probability  \t Reoccurence Interval\n')
    fprintf('------------------------------------------------------------------------------------------------\n')
    for kk = 1:n_years
        fprintf('Year %i  \t   %.4f m/s  \t   %i  \t   %f  \t   %0.1f\n',kk,annual_max(kk),rank(kk),P(kk),reoccur_interval(kk))
    end
    fprintf('------------------------------------------------------------------------------------------------\n')
    
    
    % interpolate for in between values of maximums (will return this after the plotting)
    reoccur_interval_interp = interp1(x_sorted, reoccur_interval, x_target, 'linear', 'extrap');

    if reoccur_interval_interp < 1 
        error('The reoccurence interval is less than 1 Titan year')
    else
        % Make a plot
        figure('Name','Reoccurence Interval');
        plot(reoccur_interval, x_sorted, '-k', 'LineWidth', 2); hold on;
            
        plot(reoccur_interval_interp, x_target, 'ro', 'MarkerSize', 10,'MarkerFaceColor','r', 'LineWidth', 1.5);
        xline(reoccur_interval_interp)
        yline(x_target)
    
        set(gca,'Xscale','log')
        xlabel('reoccurence interval [Titan years]');
        ylabel('surface wind speed [m/s]');
        title(sprintf('The reoccurence interval for %.2f wind speeds is %.1f Titan years\n',x_target,reoccur_interval_interp));
        grid on;
    
        
        fprintf('The reoccurence interval for %.2f wind speeds is %.1f Titan years\n',x_target,reoccur_interval_interp)
    end

    reoccur_interval = reoccur_interval_interp;

end