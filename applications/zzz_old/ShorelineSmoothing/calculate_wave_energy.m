function wave_energy = calculate_wave_energy(H,relative_angle)

    wave_energy = H.^2.*cos(relative_angle);
    if cos(relative_angle) < 0
        wave_energy = NaN;
    end
end