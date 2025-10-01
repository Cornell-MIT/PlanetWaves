
close all

[freqs,dlnf] = make_frequency_vector(Model);
dir_bins = 0:(2*pi)/Model.Dirdim:2*pi;
dir_bins(end) = [];

E_k_theta = squeeze(wn_e_spectrum{end}.E(Model.long,Model.lat,:,:));
E_k_theta = round(E_k_theta,5);
E_k_theta(E_k_theta==0) = NaN;
E_k_theta = E_k_theta';
figure;
contourf(freqs,rad2deg(dir_bins),E_k_theta)
xlabel('frequency')
ylabel('direction')
set(gca, 'XScale', 'log')
colorbar;
% caxis([0 50])
title('3 m/s Titan')
set(gca, 'XScale', 'log')
% set(gca,'ColorScale','log')
xlabel('Frequency [Hz]')
ylabel('Direction [deg]')
grid on;

%writematrix(E_k_theta,'PLANETWAVE_TITAN.txt','Delimiter','comma')