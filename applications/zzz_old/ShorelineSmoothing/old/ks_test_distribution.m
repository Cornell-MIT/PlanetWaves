function isDifferent = ks_test_distribution(data1, data2, alpha)
    % KS_TEST_DISTRIBUTION Performs Kolmogorov-Smirnov test to check if two
    % distributions are statistically different.
    %   data1, data2: Input sample data from two distributions
    %   alpha: Significance level (e.g., 0.05)
    %   isDifferent: Returns true if distributions are different, false otherwise
    
    % Perform the Kolmogorov-Smirnov test
    [h, ~] = kstest2(data1, data2, 'Alpha', alpha);
    
    % Return true if distributions are different, false otherwise
    isDifferent = logical(h);
end