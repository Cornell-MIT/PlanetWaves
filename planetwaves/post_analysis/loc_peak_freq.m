function [peak_freq,peak_freq_ind,loc_peak,dir_ind] = loc_peak_freq(energy_spectrum,gridx,gridy,model)
% location of peak frequency within energy spectrum
% INPUT: 
%   energy_spectrum = 4D energy spectrum with dimensions (Xgrid, Ygrid, Freq, Direction)
%   gridx           = location on grid in x dimension where max frequency will be determined
%   gridy           = location on grid in y dimension where max frequency will be determined
% OUTPUT:
%   peak_freq       = value of peak frequency [Hz]
%   peak_freq_ind   = index where peak frequency occurs in frequency dimension
%   loc_peak        = direction where peak frequency occurs [radians]
%   dir_ind         = index where peak frequency occurs in direction dimension

                                         
directions = 0:(2*pi)/model.Dirdim:2*pi;
directions(end) = [];

f_vec=(log(model.max_freq)-log(model.min_freq))/(model.Fdim-1);             
f = exp(log(model.min_freq)+(0:model.Fdim-1)*f_vec);                        


freq_dir = squeeze(energy_spectrum(gridx,gridy,:,:));       % 2D array of freq cols and dir rows
[~,li] = max(freq_dir(:));
[peak_freq_ind,dir_ind,] = ind2sub(size(freq_dir),li);

peak_freq = f(peak_freq_ind);
loc_peak = directions(dir_ind);
% figure;
% plot(f,tot_freq,'LineWidth',2);
% hold on;
% xline(peak_freq,'-k','LineWidth',2)
% set(gca,'XScale','log')
% xlabel('freq [hz]')
% ylabel('energy spectrum integrated over direction')
% grid on
% 
% figure;
% contourf(f,directions,freq_dir)
% hold on;
% xline(peak_freq,'--','LineWidth',2)
% yline(loc_peak,'--','LineWidth',2)
% set(gca,'XScale','log')
% xlabel('freq [hz]')
% ylabel('dir [rad]')
% view(2)
% colorbar

end