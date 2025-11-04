function [theta_centers, Epdf] = make_wave_Epdf(H, theta,rho,g)
% H = time series of all waves at a point
% theta = time series of all wave directions at a point
% rho = liquid density
% g = gravity
    
    make_plot = 1;
    num_H_bins = 20;
    num_theta_bins = 360;

    % lower edge is exclusive and the upper edge is inclusive
    % except for the very first bin, which includes both edges
    H_edges = linspace(min(H), max(H), num_H_bins+1);
    theta_edges = linspace(0, 359, num_theta_bins+1);

    E = (1/8)*rho*g.*(H.^2); % wave energy = (1/8)*rho*g*Hs^2

    Epdf = zeros(num_H_bins,num_theta_bins);
    
    % build up 2D histogram based on wave height and direction
    % point at time t has a value H(t) in direction theta(t)
    for t = 1:numel(H)
   
        theta(t) = mod(theta(t), 360);
        % Find bin indices
        h_idx = discretize(H(t), H_edges);
        theta_idx = discretize(theta(t), theta_edges);
    
        if isnan(h_idx) || isnan(theta_idx)
            error('value at %i not assigned to a bin in Epdf',t)
        end
    
        % Accumulate energy into the appropriate bin
        Epdf(h_idx, theta_idx) = Epdf(h_idx, theta_idx) + E(t);
    end

    % Epdf is 1D vector for a point on the shoreline for angles 0:360
    Epdf = sum(Epdf,1);
    % normalize probability to 1
    Epdf = Epdf./sum(Epdf(:));

    theta_centers = theta_edges(1:end-1) + diff(theta_edges)/2;

    if make_plot
        figure;
        polarplot(deg2rad(theta_centers),Epdf,'-k','LineWidth',2)
        title('energy pdf')
    end

end

