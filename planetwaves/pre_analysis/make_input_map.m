function make_input_map(Planet,Model,Wind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Makes map of bathymetry used as input for the model and includes a star at the buoy location
% of interest with a red arrow in the direction of the wind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    sz = size(Model.bathy_map);
    Xrange =  0:sz(2)-1;
    Xrange = Model.gridX.*Xrange;
    Yrange = 0:sz(1)-1;
    Yrange = Model.gridY.*Yrange;
    h = pcolor(Xrange./1000,Yrange./1000,Model.bathy_map);
    set(h,'EdgeColor','none')
    hold on;
    contour(Xrange./1000,Yrange./1000,Model.bathy_map,'-k')
    [wx,wy] = pol2cart(Wind.dir,10);
    plot(Model.long*(Model.gridX/1000),Model.lat*(Model.gridY/1000),'pentagram','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20)
    quiver(Model.long*(Model.gridX/1000),Model.lat*(Model.gridY/1000), 1*wx, 1*wy, 'r', 'MaxHeadSize', 1);
    colorbar;
    title(['input bathymetry at ',Planet.name])
    xlabel('lon [km]')
    ylabel('lat [km]')
    drawnow

end