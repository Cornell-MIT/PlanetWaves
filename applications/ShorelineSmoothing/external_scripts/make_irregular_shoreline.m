function [xcoast,ycoast,topo] = make_irregular_shoreline(Nx,Ny)

maxRetry = 100;
attempts = 0;
success = false;

% Ny = 400; % Number of grid points in y direction
% 
% Nx = 400; % Number of grid points in x direction

beta = 1.6; % Negative slope of the power spectrum. 0 = white noise, more positive values are "redder" (more variance at longer wavelengths)

variance = 1; % Variance of elevation (m^2)

rfactor = 1; % ratio of depth of closed depression to relief of the noise surface

pctwet = 30; % pctwet percent of the domain will be below sea level (z = 0)

periodic = 1; % Elevations will be periodic at the boundaries (1) or not (0, default)



% rng(0)

while attempts < maxRetry && ~success
    
    try

    % create a random component
    
    noise = RedNoise(Ny,Nx,beta,variance,periodic);
    
    relief = max(noise(:)) - min(noise(:));
    
    
    
    % taper the random component to zero at the right edge (positive x direction)
    
    X = meshgrid(1:Nx,1:Ny);
    
    taper = 1 - 1/(Nx-1)*(X-1);
    
    
    
    % create a depression with a desired depth (expressed as a multiple of the noise relief) 
    
    depth = rfactor*relief;
    
    depression = -1 * depth * Hann2D(ones(size(noise))); % this has a max of zero and a min of -depth
    
    
    
    % add the depression to the noise to create the topography.
    
    topo = depression + noise.*taper;
    
    
    
    % adjust the elevations so pctwet percent of the domain is below zero elevation
    
    % Zshift = prctile(topo(:),pctwet);
    
    Zshift = interp1(linspace(1/numel(topo(:)),1,numel(topo(:))), sort(topo(:)), pctwet/100); % poor man's percentile
    
    topo = topo - Zshift;
    
   
    
    % From taylor_order_shoreline.m
    
    conn = 8;
    
    
    
    liquid = topo < 0;
    
    
    
    land = double(~liquid);
    
    
    
    [B,L,N,A] = bwboundaries(land,conn);
    
    
    
    % As long as there isn't liquid on the boundary of the grid, the first
    
    % object should be the main landmass. The holes of that object will be the
    
    % lakes/seas, and the object children of the main object will be islands --
    
    % and those will have already had their coasts traced!
    
    
    
    k = 1; % the index of the object that is the main landmass; we assume it is 1
    
    
    
    lakes = find(A(:,k))'; % the indices of the holes that are lakes/seas
    
    
    
    % We know that each of the points on the boundary of each lake has at least
    
    % one neighbor that is land. So all we have to do is find one of those land
    
    % points, and then we can use bwtraceboundary to trace the 4-connected land
    
    % coastline.
    
    
    
    nextrow = [1,0,0; 1,0,-1; 0,0,-1];
    
    nextcol = [0,-1,-1; 0,1,0; 1,1,0];
    
    
    
    coasts = cell(size(B)); % this will hold the ordered (CW) coastlines of the lakes/seas
    
    
    
    for l = lakes 
    
        bdy = B{l};
    
        r1 = bdy(1,1);
    
        c1 = bdy(1,2);
    
        r2 = bdy(2,1);
    
        c2 = bdy(2,2);
    
        
    
        % We know that point 2 is clockwise from point 1 on the 8-connected
    
        % edge of the "water". So we go CCW one increment (45 degrees) in the
    
        % 8-connected neighborhood of point 1, and that's our land point.
    
        dr = r2-r1; % the difference in row index
    
        dc = c2-c1; % the difference in column index (note that columns increase down!)
    
        ridx = 2+dr;
    
        cidx = 2+dc;
    
        rowCCW = r1 + nextrow(ridx,cidx);
    
        colCCW = c1 + nextcol(ridx,cidx);
    
        
    
        coasts{l} = bwtraceboundary(land,[rowCCW,colCCW],'S',conn,Inf,'counterclockwise');
    
    end
    
    
    
    % Now we add in the boundaries of the islands, which we already have from
    
    % bwboundaries
    
    islands = 2:N; % the indices of the objects that are islands
    
    coasts(islands) = B(islands); % The non-empty elements of the cell array coasts 
    
                                  % should now have CW-ordered coasts for the land, not the liquid.
    
    
    
    % Find the longest lake coastline; that's probably the one we want
    
    maxlen = 0;
    
    cind = 0;
    
    for c = lakes
    
        coast = coasts{c};
    
        clength = length(coast(:,1));
    
        if clength > maxlen
    
            maxlen = clength;
    
            cind = c;
    
        end
    
    end
    
    
    
    coast = coasts{cind};
    
    xcoast = coast(:,2);
    
    ycoast = coast(:,1);
    
    
    
    % let's see what we got
    
    
    
    % figure
    
    % imagesc(L); axis image
    
    % hold on
    
    % c = cind;
    
    % % for c = [islands,lakes]
    
    %     coast = coasts{c};
    
    %     plot(xcoast,ycoast,'w','LineWidth',2);
    
    % % end
    
        
    figure;
    
    imagesc(1:Nx,1:Ny,topo)
    
    axis equal tight
    
    set(gca,'ydir','normal')
    
    colorbar
    
    hold on
    
    contour(1:Nx,1:Ny,topo,[0 0],'w')
    
    plot(xcoast,ycoast,'m','LineWidth',2);
    success = true;
    catch ME
        attempts = attempts + 1;
        warning('Shoreline did not close. Attempt %i of %i.\n %s\n',attempts,maxRetry,ME.message)
        
    end
end

end

