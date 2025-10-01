function fh = plot_wave_grid(x,y,width,height,Model,WAVES)

POI = 1;

num_colors = 100;
color1 = [229, 211, 82] / 255;
color2 = [172, 57, 49] / 255;

mycmap = [linspace(color1(1), color2(1), num_colors)', ...
                       linspace(color1(2), color2(2), num_colors)', ...
                       linspace(color1(3), color2(3), num_colors)'];


fh = figure;
plot(x,y,'-k')
hold on

for ii = 1%:Model.LonDim
    for jj = 1%Model.LatDim

        % boundaries of grid cell
        x_lower = min(x) + (ii-1) * width;
        x_upper = min(x) + ii * width;
        y_lower = min(y) + (jj-1) * height;
        y_upper = min(y) + jj * height;

        % PLOT WAVE GRID ON SUB-GRID ATTACK ANGLE
        % using relative size compared to max for coloring grid
        maxH = max(max(WAVES.H));
        H_color = round((WAVES.H(jj,ii)/maxH)*100);
        if H_color == 0
            H_color = 1;
        end
        % Check if the point of interest (POI) is in the current rectangle
            if x(POI) >= x_lower && x(POI) <= x_upper && y(POI) >= y_lower && y(POI) <= y_upper
                % Highlight the rectangle in red
                myedgecolor = 'r';
            else
                myedgecolor = 'k';
            end
            
            % Draw regular grid
            if isnan(H_color)
                rectangle('Position', [x_lower, y_lower, width, height], ...
                          'FaceColor', [1 1 1], ...
                          'EdgeColor', myedgecolor, ...
                          'FaceAlpha', 0); 
            else
                rectangle('Position', [x_lower, y_lower, width, height], ...
                          'EdgeColor', myedgecolor, ...
                          'FaceColor', mycmap(H_color, :), ...
                          'FaceAlpha', 0.5);
            end
            
    end
end

end