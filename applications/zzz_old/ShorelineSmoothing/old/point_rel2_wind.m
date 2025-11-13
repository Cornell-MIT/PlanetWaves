function my_angle = point_rel2_wind(x,y,wind_dir)
    polyin = polyshape(x,y);
    [xcent,ycent] = centroid(polyin);
    v1 = [cos(wind_dir), sin(wind_dir)];
    
    for pt = 1:numel(x)
        v2 = [(x(pt) - xcent), (y(pt) - ycent)];
        my_angle(pt) = angle_between(v1,v2);
    end



end
