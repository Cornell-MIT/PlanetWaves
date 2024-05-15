function RESULTS = analyze_buoy_data(filen)

    year = str2double(filen(end-3:end));
    buoy = filen(end-9:end-5);

    
    data_cadence = 6;
    window_size = 10;
    u_t = 3;
    dir_t = 15;
    g_t = 1.5;

    [qi, umag, gust, dir, waveht, avgu, avght, avgdir, std_u, std_ht, none_avail] = find_quiet_GREATLAKES(filen, data_cadence, window_size, u_t, dir_t, g_t);
    
    RESULTS = [avgu' avgdir' avght'];
    RESULTS = rmmissing(RESULTS);

    loc = fullfile('..','..','..','WindFetchLS_45004.csv');
    dir_fetch = csvread(loc, 1, 0);
    dir_fetch = dictionary(dir_fetch(:,1), dir_fetch(:,2));

    if ~none_avail
   
        avgfetch = dir_fetch(round(RESULTS(:,2)));
    
        RESULTS = [RESULTS avgfetch];
    
        RESULTS = table(round(RESULTS(:,1),1), round(RESULTS(:,2)), round(RESULTS(:,3),1), RESULTS(:,4), 'VariableNames', {'WSPD','WDIR','SIGHT','FETCH'});
    
        RESULTS = unique(RESULTS);
        RESULTS = [table(year.*ones(size(RESULTS(:,1))), 'VariableNames', {'YEAR'}) RESULTS];
    
        %writetable(RESULTS, sprintf('%s_%d.csv', buoy, year));
    
        

    end
end
