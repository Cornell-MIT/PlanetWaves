function [Z_value, row, col] = findNearestZ(X, Y, Z, xq, yq)
% X = Xmesh of Z
% Y = Ymesh of Z
% Z = 3D value (e.g. water depth)
% xq = x value of POI
% yq = y value of POI

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

    % If the selected cell is NaN check nearest valid neighbor and chose closest
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

    % If no valid neighbor set to NaN
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