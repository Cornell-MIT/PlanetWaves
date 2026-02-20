function RESULTS = analyze_buoy_data(filen,viable_winds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION OBJECTIVE: 
% Enforces rules for what wind data to collect 
% INPUTS:
%   filen               = Name of data file with buoy data
%   viable_winds        = criteria for relatively constant wind climate
%       .u              = magnitude of variability in wind magnitude acceptable [m/s]
%       .dir            = magnitude of variability in wind direction acceptable [deg]
%       .gust           = percent gustiness acceptable [%]
% OUTPUTS: 
%   RESULTS             = Table (Year | wind speed | wind direction | wave height ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: U.G. Schneck (schneck.una@gmail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % details about data files
    data_cadence = 6;                                                      % number of data points per hour 
    window_size = 10;                                                      % size of window of interest (controlled by the fetch and speed of the waves of interest
    year = str2double(filen(end-3:end));
    % buoy = filen(end-9:end-5);

    [~, ~, ~, ~, ~, avgu, avght, avgdir, ~, ~, none_avail] = find_quiet_GREATLAKES(filen, data_cadence, window_size, viable_winds);
    
    RESULTS = [avgu' avgdir' avght'];
    RESULTS = rmmissing(RESULTS);

    % read in associated fetch for each wind direction (from find_fetch.py) 
    loc = fullfile('..','data','Earth','WindFetchLS_45004.csv');
    dir_fetch = csvread(loc, 1, 0);
    dir_fetch = dictionary(dir_fetch(:,1), dir_fetch(:,2)); 

    if ~none_avail
   
        avgfetch = dir_fetch(round(RESULTS(:,2)));
    
        RESULTS = [RESULTS avgfetch];
        RESULTS = table(round(RESULTS(:,1),1), round(RESULTS(:,2)), round(RESULTS(:,3),1), RESULTS(:,4), 'VariableNames', {'WSPD','WDIR','SIGHT','FETCH'}); % wind speed [m/s] | wind direction [rad] | sig wave height [m] | fetch [m]
        RESULTS = unique(RESULTS);
        RESULTS = [table(year.*ones(size(RESULTS(:,1))), 'VariableNames', {'YEAR'}) RESULTS];
        %writetable(RESULTS, sprintf('%s_%d.csv', buoy, year));
    
    end
end
