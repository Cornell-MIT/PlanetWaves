function [smoothed_waveheights] = smooth_waves(waveheights)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMOOTHES SOME OF THE NUMERICAL ARTIFICATS AT HIGH WIND SPEEDS
% Useful for calculating entrainment which is sensitive to small deviations in wave height
% INPUT: 
%   waveheights          : vector of wave heights
% OUTPUT:
%   smoothed_waveheights : vector of waveheights with numerical artifacts removed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    smoothed_waveheights = smoothdata(waveheights,"rloess","SmoothingFactor",0.6);
    smoothed_waveheights(smoothed_waveheights<0) = 0;
    smoothed_waveheights(waveheights<=0) = 0;
    %ind_ceil = find(smoothed_waveheights>max(waveheights),1,'first');
    %smoothed_waveheights(ind_ceil:end) = max(waveheights);
    % dip_ind = diff(smoothed_waveheights);
    % dip_ind = find(dip_ind<0);
    % for ii = 1:numel(dip_ind)
    %     smoothed_waveheights(dip_ind(ii)) = smoothed_waveheights(dip_ind(ii)+1);
    % end

end