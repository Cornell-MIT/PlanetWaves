function make_input_map(Planet,Model,Wind)

    figure;
    imagesc(Model.bathy_map)
    hold on;
    contour(Model.bathy_map,'-k')
    [wx,wy] = pol2cart(Wind.dir,1);
    plot(Model.long,Model.lat,'pentagram','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20)
    quiver(Model.long,Model.lat, 1*wx, 1*wy, 'r', 'MaxHeadSize', 1);
    colorbar;
    title(['input bathymetry at ',Planet.name])
    drawnow

end