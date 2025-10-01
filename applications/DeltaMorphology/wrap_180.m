function new_angle = wrap_180(angle)
    
    new_angle = mod(angle + 180, 360) - 180;

end