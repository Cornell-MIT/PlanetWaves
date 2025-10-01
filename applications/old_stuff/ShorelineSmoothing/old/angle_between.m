function myangle = angle_between(vec1,vec2)
    
    % normalize the vectors
    v1 = vec1/norm(vec1);
    v2 = vec2/norm(vec2);

    % find dot product
    dot_product = dot(v1,v2);
    % find cross product
    cross_product = det([v1; v2]);

    % find angle between difference vector
    myangle = atan2(cross_product, dot_product);

    % wrap between 0 and 2pi
    if myangle < 0
        myangle = myangle + 2*pi;
    end

end
