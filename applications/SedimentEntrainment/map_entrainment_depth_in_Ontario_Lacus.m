function map_entrainment_depth_in_Ontario_Lacus(entrain_depth_m)
% Creates a map of the entrainment contours in Ontario Lacus for ArcGIS. 
% This function assumes two bathymetric slopes based on section A and L in 
% Hayes+2010 1 to 2e-3. The bathymetry is assumed to have a constant slope at all points
% where depth is given by the largest fetch for given slope.
% INPUTS
%      entrain_depth_m = depth sediment is predicted to begin to entrain [m] (e.g., sand entrains to a depth of 6 m for 1.42 m/s)
% OUTPUTS
% no variables are output but two files are created and one figure is plotted by default
%       Figure (1)    = Ontario Lacus shoreline with entrainment contours for range of slopes [deg]
%       A_basin.csv   = (Lon,Lat) contours for slope from section A in Hayes+2010 for mapping back 
%                        into ArcGIS on RADAR images
%       L_basin.csv   = (Lon,Lat) contours for slope from section L in Hayes+2010 for mapping back 
%                        into ArcGIS on RADAR images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    make_plot = 0; % make additional figures to check on conversion between 
                   % deg and meters for troubleshooting (can be left off if not debugging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2) Get mapped Ontario Lacus shoreline [deg] -> [m]
    addpath(fullfile('..','..','data/Titan/TitanLakes/shoreline'))
    load('ontariolacus_shoreline.mat','lati','lon')
    lon_shore = lon;           % mapped longitude from ArcGIS mapping [deg]
    lat_shore = lati;          % mapped longitude from ArcGIS mapping [deg]
    
    if lat_shore(1) ~= lat_shore(end)           % close the shoreline
        lon_shore = [lon_shore ;lon_shore(1)];   
        lat_shore = [lat_shore ;lat_shore(1)];
    end
    
    % deg -> m
    [x_shore, y_shore, distances,ref_pt] = titan_geodesic_distances(lat_shore,lon_shore);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (3) Create synthetic bathymetries for possible slopes
    used_saved_data = false; % false if dont want to run bathtub function to make new bathymetry each time
    
    A_slope = 2e-3; % section A from altimeter Hayes2010
    L_slope = 1e-3; % section L from altimeter Hayes2010
    
    % saved data to make this run faster but can run functions otherwise
    if used_saved_data
        addpath(fullfile('..','..','data/Titan/TitanLakes/Bathymetries/bathtub_bathy')) % contains function to make synthetic bathymetry assuming an arbitrary slope
        [A_Xmesh,A_Ymesh,A_zDep] = make_bathtub_lake(A_slope,[x_shore y_shore]);
        [L_Xmesh,L_Ymesh,L_zDep] = make_bathtub_lake(L_slope,[x_shore y_shore]);
    else
        load('past_runs/A_slope.mat','A_Xmesh','A_Ymesh','A_zDep')
        load('past_runs/L_slope.mat','L_Xmesh','L_Ymesh','L_zDep')
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (5) Get path for altimetry taking during T49 at Ontario Lacus [deg] -> [m]
    addpath(fullfile('..','..','data/Titan/TitanLakes/altimeter_pass'))
    burst_data = readtable('T49_OL.xls');
    x_path = burst_data.x;
    y_path = burst_data.y;
    x_path = unwrap_longitude(flip(x_path)); % avoids jump and needs to be flipped because unwrap uses the first phase
    x_path = flip(x_path); % flip it back baby
    
    % deg -> m
    [x_49,y_49,~] = titan_geodesic_distances(y_path,x_path,ref_pt);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (6) Plot shoreline, assumed bathymetry, and T49 altimetry pass
    if make_plot
        figure('Name','Bathymetry of Ontario Lacus with T49 burst footprint');
        subplot(2,1,1)
        plot(x_shore,y_shore,'-m','LineWidth',3)
        hold on
        contourf(A_Xmesh,A_Ymesh,A_zDep)
        plot(x_49,y_49,'-r','LineWidth',2)
        colorbar
        title(sprintf('A slope = %.1e',A_slope))
        subplot(2,1,2)
        plot(x_shore,y_shore,'-m','LineWidth',3)
        hold on
        contourf(L_Xmesh,L_Ymesh,L_zDep)
        plot(x_49,y_49,'-r','LineWidth',2)
        colorbar
        title(sprintf('L slope = %.1e',L_slope))
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (7) Plot entrainment depth in assumed bathymetry
    if make_plot
        figure('Name','Sand Entrainment in Ontario Lacus');
        plot(x_shore,y_shore,'-k','LineWidth',2)
        hold on
        contour(A_Xmesh,A_Ymesh, A_zDep,[entrain_depth_m entrain_depth_m],'--k','LineWidth',2)
        contour(L_Xmesh,L_Ymesh, L_zDep,[entrain_depth_m entrain_depth_m],'-k','LineWidth',2)
        box on;
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        plot(x_49,y_49,'-r','LineWidth',2)
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (8) Retrieve coordinates for entrained sand to map back onto RADAR image
    
    % get coordinates for sand entrainment in both basins
    [A_x, A_y] = get_contour_xy(A_Xmesh,A_Ymesh,A_zDep,entrain_depth_m);
    [L_x, L_y] = get_contour_xy(L_Xmesh,L_Ymesh,L_zDep,entrain_depth_m);
    % m -> deg
    [A_lat,A_lon] = titan_xy_to_latlon(A_x{1},A_y{1},ref_pt);
    [L_lat,L_lon] = titan_xy_to_latlon(L_x{1},L_y{1},ref_pt);
    
    figure('Name','Entrainment of Sand in Ontario in deg (for ArcGIS)')
    plot(lon_shore,lat_shore,'-r','DisplayName','Shoreline')
    hold on
    plot(A_lon,A_lat,'--k','DisplayName',sprintf('Sand entrained for %.1e slope',A_slope))
    plot(L_lon,L_lat,'-k','DisplayName',sprintf('Sand entrained for %.1e slope',L_slope))
    plot(x_path,y_path,'*b','DisplayName','T49 Altimeter Path')
    legend('show','Location','best')
    title(sprintf('Sand entrained at %.1f m',entrain_depth_m))
    grid on
    box on
    
    A_basin = [A_lon' A_lat'];
    L_basin = [L_lon' L_lat'];
    
    % write into CSV to plot in ArcGIS
    % writematrix(A_basin,'A_basin.csv','Delimiter',',')
    % writematrix(L_basin,'L_basin.csv','Delimiter',',')
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [x, y, distances, ref_point] = titan_geodesic_distances(lat_deg, lon_deg, ref_point)
    % Computes UTM-like coordinates and distances on Titan
    %   INPUT:
    %       lat_deg   = array of latitudes [degrees]
    %       lon_deg   = array of longitudes [degrees]
    %       ref_point = (lat=0,lon=0) point on map [radians]
    %   OUTPUT
    %       x         = x coordinate in meters relative to ref_point
    %       y         = y coordinate in meters relative to ref_point
    %       distances = cumulative great-circle distance in meters
    %       ref_point = (lat=0,lon=0) point on map [radians]  -- returns to be reused in future calls
    
        % Titan's mean radius (m)
        R = 2574.36 * 1000;  % Zebker+2009
    
        x = zeros(size(lat_deg));
        y = zeros(size(lat_deg));
    
        if numel(lat_deg) ~= numel(lon_deg)
            error('Latitude and longitude arrays must be the same length.');
        end
        if sum(isnan(lon_deg)) > 0 || sum(isnan(lat_deg)) > 0
            [lon_deg,lat_deg] = clean_nans(lon_deg,lat_deg);
        end
        if numel(lat_deg) < 2
            warning('Only one point specified: (lat,lon) = (%f , %f)\n',lat_deg,lon_deg)
            x = 0; y = 0; distances = 0; 
            return;
        end
    
        % Convert to radians
        lat_rad = deg2rad(lat_deg(:));
        lon_rad = deg2rad(lon_deg(:));
        
        % initialize the reference point if not specified and put in radians
        if nargin < 3 || isempty(ref_point)
            lat0 = lat_rad(1);
            lon0 = lon_rad(1);
            ref_point = [lat0, lon0];
        else
            lat0 = ref_point(1);
            lon0 = ref_point(2);
        end
    
        
        % compute distance between points
        for i = 1:length(lat_rad)
            % distane in radians for consecutive lat and lon pts
            dlat = lat_rad(i) - lat0;
            dlon = lon_rad(i) - lon0;
            % equirectangular projection (UTM-like)
            x(i) = R * dlon * cos((lat_rad(i) + lat0)/2);
            y(i) = R * dlat;
        end
    
        % Haversine segment distances for distance on sphere
        dlat = diff(lat_rad); dlon = diff(lon_rad);
        lat1 = lat_rad(1:end-1); lat2 = lat_rad(2:end);
        a = sin(dlat/2).^2 + cos(lat1).*cos(lat2).*sin(dlon/2).^2;
        c = 2 * atan2(sqrt(a), sqrt(1 - a));
        segment_distances = R * c;
    
        % Cumulative distance
        distances = [0; cumsum(segment_distances)];
    
        function [x,y] = clean_nans(x,y)
            % removes pair of coordinates if one or both values are NaN
            x(isnan(y)) = [];
            y(isnan(y)) = [];
            y(isnan(x)) = [];
            x(isnan(y)) = [];
        end
    
    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function unwrapped_lon = unwrap_longitude(lon_deg)
    % avoid jumps at ±180° for longitude
    %   Input:
    %       lon - vector of longitude values (in degrees)
    %   Output:
    %       unwrapped_lon - adjusted longitudes for smooth plotting
    
        % makes column vector
        lon_deg = lon_deg(:);
        
        % Convert to radians for unwrapping
        lon_rad = deg2rad(lon_deg);
        
        % Unwrap in radians
        unwrapped_rad = unwrap(lon_rad);
        
        % Convert back to degrees
        unwrapped_lon = rad2deg(unwrapped_rad);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [xc, yc] = get_contour_xy(X, Y, Z, level)
    % Get the (x,y) coordinates of the contours from (Xmesh,Ymesh,Depth) at some level
    %   INPUT
    %       Xmesh = X meshgrid for depth of size Z
    %       Ymesh = Y meshgrid for depth of size Z
    %       Z     = 2D array of depths
    %       level = depth contour I want to extract
    %   OUTPUT
    %       (xc,yc) = cell containing (x,y) coordinates for contour level 
    
        C = contourc(X(1,:), Y(:,1), Z, [level level]);
    
        xc = {};
        yc = {};
        k = 1;
    
        i = 1;
        while i < size(C, 2)
            this_level = C(1, i);
            n_points = C(2, i);
    
            if this_level == level
                x = C(1, i+1:i+n_points);
                y = C(2, i+1:i+n_points);
                xc{k} = x;
                yc{k} = y;
                k = k + 1;
            end
    
            i = i + n_points + 1;
        end
    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [lat_deg, lon_deg] = titan_xy_to_latlon(x, y, ref_point)
    % converts local Titan UTM in meters (x,y) pairs back to (lat,lon) in degrees
    % INPUTS
    %       x         = x coordinate in meters relative to ref_point
    %       y         = y coordinate in meters relative to ref_point
    %       ref_point = (lat=0,lon=0) point on map [radians]
    % OUTPUTS
    %       lat_deg   = array of latitudes [degrees]
    %       lon_deg   = array of longitudes [degrees]
    
        % Titan radius (m)
        R = 2574.36 * 1000;  % Zebker+2009
    
        lat0 = ref_point(1);
        lon0 = ref_point(2);
    
        % Convert back
        lat_rad = y / R + lat0;
        lon_rad = x ./ (R * cos((lat_rad + lat0)/2)) + lon0;
    
        % Convert to degrees
        lat_deg = rad2deg(lat_rad);
        lon_deg = rad2deg(lon_rad);
    end

end