function f = make_lake_wind_plot(x,y,wind_angle,my_plot_x,my_plot_y)
% angle in degrees waves are going to (CCW from East)

    f = figure;
    plot(x, y, 'k-', 'LineWidth', 2);
    hold on;
    theta_rad = deg2rad(wind_angle);
    dx = cos(theta_rad);
    dy = sin(theta_rad);
    arrow_scale = range(x)/10;
    quiver(mean(x), mean(y), dx*arrow_scale, dy*arrow_scale, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    text(mean(x)+dx*arrow_scale, mean(y)+dy*arrow_scale, 'Wind', 'Color','r','FontSize',12);
    drawnow
    



end

