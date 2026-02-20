function [u,h,angle] = table_quiet_times(wind_var,wind_dir_var,gust_var)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION OBJECTIVE: 
% make a table of all of the relative constant wind climate times in the buoy data
% This is a wrapper function for compare_model.m that opens, reads, and saves relavent data.
% INPUTS: 
%   wind_var     = magnitude of wind variability allowed [m/s]
%   wind_dir_var = magnitude of wind direction variability allowed [degree]
%   gust_var     = percent of gustiness allowable
% OUTPUTS: 
%   u            = wind speed [m/s]
%   h            = wave height [m]
%   angle        = wind direction [rad]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: U.G. Schneck (schneck.una@gmail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


years = 2002:2022; 

data = cell(1, length(years));
u = data; h = data; angle = data;

viable.u = wind_var;
viable.dir = wind_dir_var;
viable.gust = gust_var;

for i = 1:length(years)
    data{i} = analyze_buoy_data(['45004h' num2str(years(i))],viable);
end


% make a nice table of quiescent wind climates 
for i = 1:length(years)

    if ~isempty(data{i})
        u{i} = data{i}.WSPD;
        h{i} = data{i}.SIGHT;
        angle{i} = data{i}.WDIR;
    else
        u{i} = NaN;
        h{i} = NaN;
        angle{i} = NaN;
    end
end


end
