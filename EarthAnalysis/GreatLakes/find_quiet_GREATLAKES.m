function [quiet_time,umag,gust,dir,waveht,avgu,avght,std_u,std_ht] = find_quiet_GREATLAKES(filename,data_cadence,window_size,u_threshold,direction_threshold,gust_threshold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION OBJECTIVE: 
% Finds a period of quiet where the wind speed is relatively constant, the gusts are relatively similar to the wind speed,
% and the direction the wind does not vary hugely.
% INPUTS:
%   filename            = Name of .txt file with Great Lake Buoy Table (e.g. for Bouy 45004 find it here: https://www.ndbc.noaa.gov/station_history.php?station=45004)
%   data_cadence        = number of data points per hour 
%   window_size         = size of window of interest (controlled by the fetch and speed of the waves of interest)
%   u_threshold         = maximum allowable change in |u| within the window of time
%   direction_threshold = maximum allowable change in the direction of the wind within the window of time
%   gust_threshold      = maximum allowable fraction of gusts relative to average wind speed (e.g. 1.5 = 150% threshold)
% OUTPUTS: 
%   quiet_time          = index of a period of quiet within the buoy data that matches threshold criteria
%   umag                = cleaned time series of wind speed [m/s]
%   gust                = cleaned time series of gusts [m/s]
%   dir                 = cleaned time series of wind direction [rad]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: U.G. Schneck (schneck.una@gmail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    yr = filename(end-3:end);

    thd_dir = direction_threshold; % threshold for direction (less than 10 degree change)
    thd_g = gust_threshold; % threshold for gust (less than 150%)
    
    % LOAD IN DATA
    warning('off','MATLAB:table:ModifiedAndSavedVarnames') % suppress annoying message about changing the year name column that I don't care about
    lakedatatable = readtable(strcat(filename,'.txt'),MissingRule="omitrow",ReadVariableNames=true);
    
    window_size = data_cadence*window_size; % window size in terms of number of measurements
    
    % extract relevant data from table
    umag = lakedatatable.WSPD;
    sig_h = lakedatatable.WVHT;
    gust = lakedatatable.GST;
    if any("WDIR" == string(lakedatatable.Properties.VariableNames))
        dir = lakedatatable.WDIR;
    else
        dir = lakedatatable.WD;
    end
    waveht = lakedatatable.WVHT;
       
    % CLEAN DATA
    % remove all points where there are no wind speed measurements
    sig_h(umag==99) = NaN;
    gust(umag==99) = NaN;
    dir(umag==99) = NaN;
    umag(umag==99) = NaN;
    
    % turn all missing values (99) into NaN 
    dir(dir==999) = NaN;
    sig_h(sig_h==99) = NaN;
    umag(sig_h==99) = NaN;
    waveht(waveht==99) = NaN;
    
    % initalize vectors
    quiet_time = NaN(1,numel(umag)-window_size);
    avgu = NaN(1,numel(umag)-window_size);
    avght = NaN(1,numel(umag)-window_size);
    std_u = NaN(1,numel(umag)-window_size);
    std_ht = NaN(1,numel(umag)-window_size);

    for i = 1:length(umag) - window_size
    
       tm = i:i+window_size;
       frac_g = max(gust(tm))/mean(umag(tm),"omitmissing");
       chg_u = max(umag(tm)) - min(umag(tm));
       chg_dir = (max(dir(tm)) - min(dir(tm)));

       if chg_u <= u_threshold && chg_dir <= thd_dir && frac_g <= thd_g
           quiet_time(i) = 1;
           avgu(i) = mean(umag(tm),"omitmissing");
           avght(i) = mean(sig_h(tm),"omitmissing");
           std_u(i) = std(umag(tm),"omitmissing");
           std_ht(i) = std(sig_h(tm),"omitmissing");
       else
           quiet_time(i) = 0;
           avgu(i) = NaN;
           avght(i) = NaN;
           std_u(i) = NaN;
           std_ht(i) = NaN;
       end

    end

    if sum(isnan(quiet_time))
        disp('nan still present')
        %quiet_time = rmmissing(quiet_time);
    end

    if ~sum(quiet_time)
        mes = sprintf('%s : no data fits this criteria, change the thresholds for |u|, direction, and/or gusts \n',yr);
        fprintf(mes);
    else
        mes =  sprintf('%s : SUCCESS \n',yr);
        fprintf(mes);
    end
       


end