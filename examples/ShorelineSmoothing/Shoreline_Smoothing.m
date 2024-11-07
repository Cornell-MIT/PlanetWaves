function [Qs,psi] = shoreline_smoothing(x,y,sig_wave,theta,WIDTH,HEIGHT,Model)

Qs = NaN(size(x));     % Sediment transport rate
psi = NaN(size(x));
% Calculate Qs and diffusivity for each (x, y) point based on the same theta ordering
for i = 1:numel(x)

    % skip if POI is NaN
    if isnan(x(i))
        continue
    end
    
    % find grid containing POI (x(i), y(i)) 
    cell_x = floor((x(i) - min(x)) / WIDTH) + 1;
    cell_y = floor((y(i) - min(y)) / HEIGHT) + 1;
    
    % skip points outside coordinate grid
    if cell_x < 1 || cell_x > Model.LonDim || cell_y < 1 || cell_y > Model.LatDim
        %continue  
        if cell_y > Model.LatDim
            cell_y = Model.LatDim;
        elseif cell_x > Model.LonDim
            cell_x = Model.LonDim;
        else
            continue
        end
    end
    
    % get H, T of current cell
    H_cell = sig_wave.H(cell_y, cell_x);
    T_cell = sig_wave.T(cell_y, cell_x);
    
    if isnan(H_cell)
        [H_cell, T_cell] = check_neighbor(sig_wave, cell_x, cell_y, Model.LonDim, Model.LatDim, x(i), y(i), min(x), min(y), WIDTH, HEIGHT);
    end

    if isnan(H_cell)
        fprintf('no valid neighbors found %i\n',i)
    end
    % Calculate sediment transport rate Qs based on theta(i)
    Qs(i) = CERC(H_cell, T_cell, theta(i));
    psi(i) = shoreline_stability(H_cell, T_cell, theta(i),Model.bathy_map(cell_x,cell_y));
end

Qs(isnan(Qs)) = 0;
psi(isnan(psi)) = 0;

end