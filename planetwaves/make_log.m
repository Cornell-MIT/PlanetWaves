function TitanResults = make_log(planet,model,wind,uniflow,Etc)
% make a log and save file location for result outputs

    
    % -- create output directory for results 
    ResultsParent = 'results';
    TitanResults = strcat('wind_speed_',num2str(wind.speed)); 
    TitanResults = fullfile(ResultsParent,TitanResults);

    
    if exist(ResultsParent,'dir') ~= 7 && Etc.savedata
        disp(['Creating folder ' ResultsParent ' to store model output'])
        mkdir(ResultsParent);
    end
    if exist(TitanResults, 'dir') == 7  && Etc.savedata                        
       oldmatfiles = fullfile(TitanResults, '*.mat');                          % empties output directory from previous runs
       oldmatloc = dir(oldmatfiles);
       disp(['Previous runs already exist in ' TitanResults])
       for kk = 1:length(oldmatloc)
           basemat = oldmatloc(kk).name;
           fullmat = fullfile(TitanResults,basemat);
           fprintf(1,'Deleting previous .mat files %s\n',fullmat);
           delete(fullmat);
       end
       oldlogfile = fullfile(TitanResults,'*.txt');
       oldlogloc = dir(oldlogfile);
       for kk = 1:length(oldlogloc)
           basetxt = oldlogloc(kk).name;
           fulltxt = fullfile(TitanResults,basetxt);
           fprintf(1,'Deleting previous log files %s\n',fulltxt);
           delete(fulltxt);
       end
    elseif Etc.savedata
        disp(['Creating subdirectories in ' ResultsParent ' to store model output'])
	    mkdir(TitanResults)
    end
    
        % -- prepare log file for commands 
    dfile=strcat(string(datetime('now','TimeZone','local','Format','ddMMyy_HHmmss')),'_wind_speed_',num2str(wind.speed),'_RunLog.txt');
    diary(fullfile(TitanResults,dfile));
    RAII.diary = onCleanup(@() diary('off'));                                  % auto-closes logging function on error

% adds run details to console and saves to log file 
    disp('================================================================')
    disp(['Directional Wave Spectrum -- last updated: ' dir('makeWaves.m').date])
    disp(['Wind Speed(s) to Run: ' regexprep(mat2str(wind.speed),{'\[', '\]', '\s+'}, {'', '', ','}) ' m/s']);
    disp('================================================================')

    fprintf('Planet:\t\t\t\t %s\n',planet.name)
    fprintf('\trho_liquid:\t\t %.3f\n',planet.rho_liquid)
    fprintf('\tnu_liquid:\t\t %.3e\n',planet.nu_liquid)
    fprintf('\tnua:\t\t\t %.3e\n',planet.nua)
    fprintf('\tgravity:\t\t %.3f\n',planet.gravity)
    fprintf('\tsurface_temp:\t\t %.3f\n',planet.surface_temp)
    fprintf('\tsurface_press:\t\t %.3f\n',planet.surface_press)
    fprintf('\tsurface_tension:\t %.3f\n',planet.surface_tension)
    
    fprintf('\n')
    
    fprintf('Model\n')
    fprintf('\tm:\t\t\t %i\n',model.LonDim)
    fprintf('\tn:\t\t\t %i\n',model.LatDim)
    fprintf('\to:\t\t\t %i\n',model.Fdim)
    fprintf('\tp:\t\t\t %i\n',model.Dirdim)
    fprintf('\tlong:\t\t\t %i\n',model.long)
    fprintf('\tlat:\t\t\t %i\n',model.lat)
    fprintf('\tgridX:\t\t\t %i\n',model.gridX)
    fprintf('\tgridY:\t\t\t %i\n',model.gridY)
    fprintf('\tmindelt:\t\t %.2e\n',model.mindelt)
    fprintf('\tmaxdelt:\t\t %.2f\n',model.maxdelt)
    fprintf('\ttime_step:\t\t %i\n',model.time_step)
    fprintf('\tnum_time_steps:\t\t %i\n',model.num_time_steps)
    fprintf('\ttolH:\t\t\t %.2f\n',model.tolH)
    fprintf('\tcutoff_freq:\t\t %i\n',model.cutoff_freq)
    fprintf('\tmin_freq:\t\t %.2f\n',model.min_freq)
    fprintf('\tmax_freq:\t\t %.2f\n',model.max_freq)
    fprintf('\tz_data:\t\t\t %.2f\n',model.z_data)
    fprintf('\ttune_A1:\t\t %.2f\n',model.tune_A1)
    fprintf('\ttune_mss_fac:\t\t %.2f\n',model.tune_mss_fac)
    fprintf('\ttune_Sdt_fac:\t\t %.3f\n',model.tune_Sdt_fac)
    fprintf('\ttune_Sbf_fac:\t\t %.3f\n',model.tune_Sbf_fac)
    fprintf('\ttune_cotharg:\t\t %.2f\n',model.tune_cotharg)
    fprintf('\ttune_n:\t\t\t %.2f\n',model.tune_n)
    fprintf('\tbathy_map (deepest):\t %.2f\n',max(max(model.bathy_map)))
    
    allSame = all(all(model.bathy_map == model.bathy_map(1,1)) == 1);
    if allSame
        slopex = 0;
        slopey = 0;
    else
        slopex = max(max(model.bathy_map))/(model.LonDim*model.gridX);
        slopey = max(max(model.bathy_map))/(model.LatDim*model.gridY);
    end
    
    
    fprintf('\tbathy_map (slopex):\t %.2e\n',slopex)
    fprintf('\tbathy_map (slopey):\t %.2e\n',slopey)
    
    
    fprintf('\n')
    
    fprintf('Wind\n')
    [fprintf('\tspeed: \t\t\t'),fprintf(' %.2f,', wind.speed(1:end-1)), fprintf(' %.2f\n', wind.speed(end))];
    fprintf('\tdir: \t\t\t %.2f',wind.dir)
    
    fprintf('\n')
    
    fprintf('Uniflow\n')
    fprintf('\tEast:\t\t\t %.2f\n',uniflow.East)
    fprintf('\tNorth:\t\t\t %.2f\n',uniflow.North)
    
    fprintf('\n')
    
    fprintf('Etc\n')
    fprintf('\tshowplots:\t\t %i\n',Etc.showplots)
    fprintf('\tsavedata:\t\t %i\n',Etc.savedata)

disp('================================================================')

end