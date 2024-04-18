clc
clear
close all
load('LS_uniform_depth_45004_04162024')

for k = 1:numel(test_speeds)
    buoy_waves(k) = htgrid{k}{end}(Model.long,Model.lat);
end

figure('units','normalized','outerposition',[0 0 1 1])
for speed = 1:numel(test_speeds)
    
    subplot(1,3,3)
    plot(test_speeds,buoy_waves,'-sb','LineWidth',1,'MarkerFaceColor','b')
    hold on
    plot(test_speeds(speed),htgrid{speed}{end}(Model.long,Model.lat),'-sr','LineWidth',1,'MarkerFaceColor','r')
    hold off;
    xlabel('u [m/s]')
    ylabel('H_{sig} [m]')
    grid on;

    ax1 = subplot(1,3,[1,2]);
    imagesc(1:Model.m,1:Model.n,htgrid{speed}{end}')
    colormap linspecer
    ylabel('longitude [km]')
    xlabel('latitude [km]')
    title(sprintf('u = %i m/s',test_speeds(speed)))
    c1 = colorbar;
    c1.Label.String = 'H_{sig} [m]';
    clim([0 8])
    hold on;
    
    contour(htgrid{speed}{end}','-k','LineWidth',2)

    grid on
    new_xtick = get(gca, 'XTick')*(Model.gridX)/1000;
    new_ytick = get(gca, 'YTick')*(Model.gridY)/1000;
    set(gca, 'XTick',  get(gca, 'XTick'), 'XTickLabel', arrayfun(@(x) sprintf('%d', x), new_xtick, 'UniformOutput', false));
    set(gca, 'YTick',  get(gca, 'YTick'), 'YTickLabel', arrayfun(@(y) sprintf('%d', y), new_ytick, 'UniformOutput', false));
    
    
    [wx,wy] = pol2cart(Wind.dir,1);
    plot(Model.lat,Model.long,'or','MarkerFaceColor','r')

    for i = 1.5:2:10
        for j = 1:2:10
            quiver(i, j, (speed/5)*wx, (speed/5)*wy, 'k', 'MaxHeadSize', 1);
        end
    end
    set(ax1,'Ydir','normal')
    set(ax1,'Xdir','reverse')
    drawnow
    if speed == 1
        gif('uniform_depth_LS_45004.gif','DelayTime',1)
    else
        gif
    end
end

