function [theta_centers, Epdf] = make_wave_Epdf(H, theta,T)
% INPUT
% H = time series of all waves heights
% T = time series of all waves periods
% theta = time series of all wave directions 
% OUTPUT
% Epdf = wave energy matrix 
% theta_centers = wave directions associated with wave heights in Epdf 

    num_H_bins = 20;
    num_theta_bins = 360;

    H(isnan(H)) = 0;
    % lower edge is exclusive and the upper edge is inclusive
    % except for the very first bin, which includes both edges
    H_edges = linspace(min(H), max(H), num_H_bins+1);
    theta_edges = linspace(0, 360, num_theta_bins+1); % e.g., bin 1 = 0 to 1 degree, bin 360 = 359 to 360 deg

    E = H.^(12/5).*T.^(1/5); % equation 3 Jaap supplement
    Epdf = zeros(num_H_bins,num_theta_bins);
    
    % build up 2D histogram based on wave height and direction
    % point at time t has a value H(t) in direction theta(t)
    for t = 1:numel(H)
   
        % Find bin indices
        h_idx = discretize(H(t), H_edges);
        theta_idx = discretize(theta(t), theta_edges);
    
        if isnan(h_idx) || isnan(theta_idx)
            error('value at %i not assigned to a bin in Epdf',t)
        end
    
        % Accumulate energy into the appropriate bin
        Epdf(h_idx, theta_idx) = Epdf(h_idx, theta_idx) + E(t); % [wave heights x wave directions]
    end

    % Epdf is 1D vector for a point on the shoreline for angles 0:360
    Epdf = sum(Epdf,1,'omitmissing'); % sum over H index so just a function of theta
    Epdf(isnan(Epdf)) = 0;
    Epdf = Epdf./sum(Epdf); % normalize PDF to sum to one
   

    theta_centers = theta_edges(1:end-1);  % bin centers: 0,1,2,..,359

    theta_centers = wrapTo180(theta_centers);

end

