function [binned_wave_energy,binned_direction] = make_wind_rose(wind_speed,wind_direction,bin_size)
% 
% INPUT
% wind_speed : magnitude of wind
% wind_direction : angle (0 = West -> East) and increasing CCW
% OUTPUT
% binned_wave_energy : PDF of wave energy from a direction 
% binned_direction : associated direction for PDF of binned_wave_energy (0 = West -> East, increasing CCW)

    make_plot = 1;

    wave_energy = wind_speed.^2;

    bin_edge = 0:bin_size:360;
    binned_direction = bin_edge(1:end-1) + bin_size/2;
    nbins = numel(binned_direction);

    binned_wave_energy = zeros(1,nbins);

    wind_direction = mod(wind_direction,360);

    for i = 1:nbins
        bin_idx = wind_direction >= bin_edge(i) & wind_direction < bin_edge(i+1);
        %binned_wave_energy(i) = numel(wave_energy(bin_idx))/numel(wave_energy);
        binned_wave_energy(i) = sum(wave_energy(bin_idx),'all','omitmissing')/sum(wave_energy,'all','omitmissing');
    end


    if make_plot
        theta = deg2rad(binned_direction);
        figure;
        polarplot([theta, theta(1)], [binned_wave_energy, binned_wave_energy(1)], '-k', 'LineWidth', 2);
        rlim([0, max(binned_wave_energy)*1.1])
        thetalim([0 360])
        ax = gca;
        ax.ThetaZeroLocation = 'right';
        ax.ThetaDir = 'counterclockwise';
        title('Wave Energy')
    end


end