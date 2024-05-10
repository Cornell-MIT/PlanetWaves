function ht_wind = plot_spectrogram(file_loc,which_time)

% plot_spectrogram(fullfile('LakeSuperior_Deepwater','wind_speeds'),10);
addpath(file_loc)

mytime = strcat('_',string(which_time),'.mat');
which_time = strcat('*_',string(which_time),'.mat');
spectra_files = dir(fullfile(file_loc,which_time));
run_details = dir(fullfile(file_loc,'New_Reference*'));

load(run_details.name,'freqs','p')
dir_bins = 0:(2*pi)/p:2*pi;
dir_bins(end) = [];

labelled_spectrum = containers.Map;

time_labels = strings(1,numel(spectra_files));

for i = 1:numel(spectra_files)
    load(spectra_files(i).name,'E')
    label = erase(spectra_files(i).name,'.mat');
    time_labels(i) = erase(spectra_files(i).name,mytime);
    time_labels(i) = erase(time_labels(i),'New_');
    labelled_spectrum(time_labels(i)) = E;
    aa = labelled_spectrum(time_labels(i));
    sz_aa = size(aa);
    pvs = squeeze(aa(round(sz_aa(1)/2),round(sz_aa(2)/2),:,:));
    pvs = round(pvs,5);
    pvs(pvs==0) = NaN;
    power_vs_spread{i} = pvs;

    pp = power_vs_spread{i}';
    % figure;
    % contour(freqs,rad2deg(dir_bins),power_vs_spread{i}')
    % colormap(linspecer)
    % colorbar;
    % caxis([0 50])
    % title(wind_speeds(i))
    % set(gca, 'XScale', 'log')
    % set(gca,'ColorScale','log')
    % xlabel('Frequency [Hz]')
    % ylabel('Direction [deg]')
    % grid on;
    save(strcat(label,"_results.mat"),'freqs','dir_bins','pp');

    load(spectra_files(i).name,'ht')    
    ht_wind(i) = ht(5,2);

    % if i == 1
    %     gif('power.gif')
    % else
    %     gif;
    % end
end
end