function LST = CERC(H,T,relative_angle)
    
    for a = 1:numel(relative_angle)
        if relative_angle(a) >= -pi/2 &&  relative_angle(a) <= pi/2
            LST(a) = (H^(12/5))*(T^(1/5))*(cos(relative_angle(a))^(6/5))*sin(relative_angle(a));
        else
            LST(a) = NaN;
        end
    end
end

