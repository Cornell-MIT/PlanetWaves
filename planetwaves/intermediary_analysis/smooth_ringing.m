function smooth_A = smooth_ringing(A,H)
    
    fraction = 10;
    T = size(A,1);

    startIdx = T - floor(T / fraction) + 1;
    
    if sum(diff(H)>0)
        smooth_A = squeeze(mean(A(startIdx:T,:, :), 1));
    else
        smooth_A = squeeze(A(end,:,:));
    end

end