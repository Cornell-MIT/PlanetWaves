% clc
% clear
close all
% 
% % % HOUSEKEEPING
% addpath(fullfile('..','..','planetwaves'))  
% addpath(fullfile('..','..','planetwaves','pre_analysis'))
% addpath(fullfile('..','..','planetwaves','post_analysis'))
% addpath(fullfile('external_scripts'))  
% jet_wrap = vertcat(jet,flipud(jet));
% 
% % make a rednoise shoreline
% % num_gridX = 400;
% % num_gridY = 400;
% % [x,y,Depth] = make_irregular_shoreline(num_gridX,num_gridY);
% load('asylake1.mat')
% % make a bathtub bathymetry for it
% load('asylake1_bathtub.mat')
% 
% buoy_loc = [200 200];
% grid_resolution = [300 300];
% planet_to_run = 'Titan-OntarioLacus';
% time_to_run = 60;
% wind_direction = pi;
% 
% zDep = zDep.*1000;
% 
% zDep(:,1:400) = [];
% zDep(:,180:600) = [];
% zDep(1:440,:) = [];
% zDep(120:560,:) = [];
% 
% Xmesh(:,1:400) = [];
% Xmesh(:,180:600) = [];
% Xmesh(1:440,:) = [];
% Xmesh(120:560,:) = [];
% 
% Ymesh(:,1:400) = [];
% Ymesh(:,180:600) = [];
% Ymesh(1:440,:) = [];
% Ymesh(120:560,:) = [];
% 
% 
% figure;
% surf(Xmesh,Ymesh,zDep)
% view(2)
% colorbar
% hold on
% plot(x,y,'-r')
% plot3(buoy_loc(1),buoy_loc(2),200,'om','MarkerFaceColor','m')
% resizeFactor = 1/10;
% 
% [zDep, buoy_loc, grid_resolution, Xmesh, Ymesh, ~, ~] = degrade_depth_mesh(zDep, buoy_loc, grid_resolution, resizeFactor, Xmesh, Ymesh, x, y);
% 
% zDep = round(zDep);
% zDep(zDep==0) = NaN;
% 
% buoy_loc = [10 5];
% 
% figure;
% surf(Xmesh,Ymesh,zDep)
% view(2)
% colorbar
% hold on
% plot3(x,y,200.*ones(size(x)),'-r') 
% 
% 
% [Planet,Model,Wind,Uniflow,Etc] = initalize_model(planet_to_run,time_to_run,wind_direction,zDep,buoy_loc);
% 
% Model.gridX = grid_resolution(1);
% Model.gridY = grid_resolution(2);
% make_input_map(Planet,Model,Wind)
% 
% 
% 
% 
% Wind.speed = 1;
% Wind.dir = wind_direction;
% 
% Model = calc_cutoff_freq(Planet,Model,Wind);
% 
% [~,~,~,~,~,~,PeakWave] = makeWaves(Planet, Model, Wind, Uniflow, Etc);  
% sig_wave = invert_attribute(PeakWave);
% 
% figure;
% surf(Xmesh,Ymesh,sig_wave.H)
% view(2)
% colorbar
% hold on
% plot3(x,y,200.*ones(size(x)),'-r') 

WIDTH = (max(x)-min(x))/Model.LonDim;                                      % Width of each grid cell (# cells * width of 1 cell = width of all cells)
HEIGHT = (max(y)-min(y))/Model.LatDim;                                     % Height of each grid cell


for POI = 1:numel(x)
    
    % Find nearest Z value
    [mywave(POI), row, col] = findNearestZ(Xmesh, Ymesh, depth, x(POI), y(POI));
    title(sprintf('Depth is %f', mywave(POI)));
    % Create surface plot
    s = surf(Xmesh, Ymesh, depth, 'FaceColor', 'flat');
    view(2);
    hold on;

    alphaData = ones(size(depth)) * 0.5; 
    alphaData(row, col) = 1; % leep selected grid cell fully visible
    set(s,'FaceAlpha','flat','AlphaDataMapping','none','AlphaData',alphaData)

     % Plot trajectory and selected point
    plot3(x, y, 100 .* ones(size(x)), '-r');
    plot3(x(POI), y(POI), 100, 'om');

    % Highlight selected grid cell
    plot3(Xmesh(row, col), Ymesh(row, col), 100, 'or');

    drawnow;
    hold off;

    % Save as GIF
    if POI == 1
        gif('waves.gif');
    else
        gif;
    end
    
end

% sinu = calc_sinuosity(x,y,100);
% 
% figure
% scatter(x,y,100,sinu,'filled')
% colorbar
% title('sinuosity')
% 
% shore_angle = point_rel2_wind(x,y,wind_dir);
% 
% 
% figure
% scatter(x,y,50,rad2deg(shore_angle),'filled')
% colormap(jet_wrap)
% colorbar
% hold on
% [u_wind,v_wind] = pol2cart(wind_dir, 50);
% quiver(200,200,u_wind,v_wind)
% title('angle relative to wind')
% 
% for i = 1:numel(x)
%     dn_dt(i) = shoreline_stability(1,1,shore_angle(i),1);
% end
% 
% along_wind = find_alongwind(x,shore_angle);
% 
% [r_along,s_along] = r_squared(dn_dt(along_wind),sinu(along_wind),'along wind diffusvity','along wind sinuosity');
% [r_notalong,s_notalong] = r_squared(dn_dt(~along_wind),sinu(~along_wind),'not along wind diffusvity','not along wind sinuosity');
% 
% figure;
% h1 = histogram(sinu(along_wind),[0:0.1:1]);
% hold on
% h2 = histogram(sinu(~along_wind),[0:0.1:1]);
% legend('along wind','oppose wind')
% xlabel('sinuosity')
% 
% isDifferent = ks_test_distribution(sinu(along_wind),sinu(~along_wind),0.05);
% if isDifferent
%     title('shoreline is asym in sinuosity wrt wind (ks test 95%)')
% else
%     title('shoreline is NOT asym in sinuosity wrt wind (ks test 95%)')
% end



function [Z_value, row, col] = findNearestZ(X, Y, Z, xq, yq)
    
    if ~isequal(size(X), size(Y), size(Z))
        error('X, Y, and Z must have the same dimensions.');
    end
    
    %find closest grid cell to query point
    col = find(X(1,:) <= xq, 1, 'last');
    row = find(Y(:,1) <= yq, 1, 'last');
    
    % If the selected grid cell is not NAN
    if ~isnan(Z(row, col))
        Z_value = Z(row, col);
        return;
    end
    
    % If the selected cell is NaN, check nearest valid neighbor and chose closest
    neighbors = [-1 -1; -1 0; -1 1; 
                  0 -1;  0 1; 
                  1 -1;  1 0;  1 1];

    minDist = inf;
    newRow = row;
    newCol = col;
    foundValid = false;

    % loop through eight neighbors
    for k = 1:size(neighbors, 1)
        r = row + neighbors(k,1);
        c = col + neighbors(k,2);


        if r >= 1 && r <= size(Z,1) && c >= 1 && c <= size(Z,2)
            if ~isnan(Z(r, c)) % ignore is NaN
                % find distance to cell's center
                centerDist = (X(r, c) - xq)^2 + (Y(r, c) - yq)^2;
                if centerDist < minDist
                    minDist = centerDist;
                    newRow = r;
                    newCol = c;
                    foundValid = true;
                end
            end
        end
    end

    % If no valid neighbor was found, set to NaN
    if ~foundValid
        warning('No suitable nearest grid cells found. Assigning Z_value to NaN.');
        Z_value = NaN;
    else
        % Use the closest valid neighbor
        row = newRow;
        col = newCol;
        Z_value = Z(row, col);
    end
end
