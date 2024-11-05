function smooth_A = smooth_ringing(A,H)
% FUNCTION TO SMOOTH OUT NUMERICAL RINGING
% WILL CHECK IF THERE IS RINGING AND THEN 
% ASSIGN THE VALUE TO THE AVERAGE OF THE 
% LAST TENTH OF THE TIME SERIES
% INPUT: 
%   A        : array to check (time x M x N)
%   H        : time series of wave heights to check for ringing
% OUTPUT:
%   smooth_A : final value of A either the last time step if no ringing or the average of the last fraction of time steps

    fraction = 10;
    T = size(A,1);

    startIdx = T - floor(T / fraction) + 1;
    
    if sum(diff(H)>0)
        smooth_A = squeeze(mean(A(startIdx:T,:, :), 1));
    else
        smooth_A = squeeze(A(end,:,:));
    end

end