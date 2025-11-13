function  wave_energy = calc_wave_energy(x,y,sig_wave,WIDTH,HEIGHT,Model,Planet,relative_angle)

testing = 0;
POI = 1; % for trouble-shooting

wave_energy = NaN(size(x));

% Calculate wave energy for each (x, y) point based on the same theta ordering
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
    H_cell = sig_wave.H(cell_x, cell_y);
    T_cell = sig_wave.T(cell_x, cell_y);

    depth = Model.bathy_map(cell_x,cell_y);
    
    
    if isnan(H_cell)
        [H_cell, T_cell,depth] = check_neighbor(sig_wave, cell_x, cell_y, Model.LonDim, Model.LatDim, x(i), y(i), min(x), min(y), WIDTH, HEIGHT,Model,x,y);
    end
    
    wave_energy(i) = (1/8)*(Planet.rho_liquid)*(Planet.gravity)*(H_cell^2);
    wave_energy(i) = wave_energy(i)*cos(relative_angle(i));

    if i == POI && testing == 1 % for trouble-shooting
        disp('===========')
        fprintf('xcell: %i\n',cell_x)
        fprintf('ycell %i\n',cell_y)
        fprintf('D: %f\n',depth)
        fprintf('H: %f\n', H_cell)
        fprintf('T: %f\n',T_cell)
    end



end

wave_energy(isnan(wave_energy)) = 0;


end