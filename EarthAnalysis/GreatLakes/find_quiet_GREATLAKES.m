function [quiet_index,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filename,data_cadence,window_size,direction_threshold,gust_threshold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION OBJECTIVE: 
% Finds a period of quiet where the wind speed is relatively constant, the gusts are relatively similar to the wind speed,
% and the direction the wind does not vary hugely.
% INPUTS:
%   filename            = Name of .txt file with Great Lake Buoy Table 
%                         (e.g. for Bouy 45004 find it here: https://www.ndbc.noaa.gov/station_history.php?station=45004)
%   data_cadence        = number of data points per hour 
%   window_size         = size of window of interest (controlled by the fetch and speed of the waves of interest)
%   direction_threshold = maximum allowable change in the direction of the wind within a window of time
%   gust_threshold      = maximum allowable fraction of gusts relative to average wind speed (e.g. 1.5 = 150% threshold)
% OUTPUTS: 
%   quiet_index         = index of a period of quiet within the buoy data that matches threshold criteria
%   umag                = cleaned time series of wind speed [m/s]
%   gust                = cleaned time series of gusts [m/s]
%   dir                 = cleaned time series of wind direction [rad]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: U.G. Schneck (schneck.una@gmail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    thd_dir = direction_threshold; % threshold for direction (less than 10 degree change)
    thd_g = gust_threshold; % threshold for gust (less than 150%)
    
    % LOAD IN DATA
    warning('off','MATLAB:table:ModifiedAndSavedVarnames') % suppress annoying message about changing the year name column that I don't care about
    lakedatatable = readtable(strcat(filename,'.txt'));
    
    window_size = data_cadence*window_size; % window size in terms of number of measurements
    
    % extract relevant data from table
    umag = lakedatatable.WSPD;
    sig_h = lakedatatable.WVHT;
    gust = lakedatatable.GST;
    dir = lakedatatable.WDIR;
    waveht = lakedatatable.WVHT;
       
    % CLEAN DATA
    % remove all points where there are no wind speed measurements
    sig_h(umag==99) = NaN;
    gust(umag==99) = NaN;
    dir(umag==99) = NaN;
    umag(umag==99) = NaN;
    
    % remove all points where the wind speed is below the threshold for wave development
    % sig_h(umag<1.6) = NaN;
    % gust(umag<1.6) = NaN;
    % dir(umag<1.6) = NaN;
    % umag(umag<1.6) = NaN;
    % turn all missing values (99) into NaN 
    dir(dir==999) = NaN;
    sig_h(sig_h==99) = NaN;
    umag(sig_h==99) = NaN;
    waveht(waveht==99) = NaN;
    
    var_u = NaN(1,numel(umag)-window_size);
    avg_u = NaN(1,numel(umag)-window_size);
    var_gust = NaN(1,numel(umag)-window_size);
    avg_gust = NaN(1,numel(umag)-window_size);
    var_dir = NaN(1,numel(umag)-window_size);
    avg_dir = NaN(1,numel(umag)-window_size);
    
    for i = 1:numel(umag)-window_size
        var_u(i) = var(umag(i:i+window_size),"omitnan");
        avg_u(i) = mean(umag(i:i+window_size),"omitmissing");
        
        var_gust(i) = var(gust(i:i+window_size),"omitnan");
        avg_gust(i) = mean(gust(i:i+window_size),"omitmissing");
    
        var_dir(i) = var(dir(i:i+window_size),"omitnan");
        avg_dir(i) = mean(dir(i:i+window_size),"omitmissing");
    end
    
    
    [usortv,usorti] = sort(var_u);
    

    for i = 1:length(usortv)
        
        tm = usorti(i):usorti(i)+window_size;
        low_cand = umag(tm); % candidate for quiet time period, first use a window of time with low wind speed variability
                
        frac_g = max(gust(tm))/mean(low_cand,"omitmissing"); % max gust / avg of wind speed in time of interest
    
        if (max(dir(tm)) - min(dir(tm))) <= thd_dir && frac_g <= thd_g
            quiet_index = tm;
            break
        end
        
    
    
    end

    avgu = mean(umag(quiet_index),"omitmissing");
    avght = mean(waveht(quiet_index),"omitmissing");
    std_u = std(umag(quiet_index),"omitmissing");
    std_ht = std(waveht(quiet_index),"omitmissing");


end