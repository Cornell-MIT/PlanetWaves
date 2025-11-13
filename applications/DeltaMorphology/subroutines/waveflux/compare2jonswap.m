function compare2jonswap(wind_speed,fetch,g,Hs_model,Tp_model)
% wind_speed is [1 x num_wind], fetch, Hsmodel, Tp_model is [num_wind x num_points]

    make_plot = 1;

    wind_speed = wind_speed(:); % [num_wind x 1] for broadcasting

    % JONSWAP approximation
    Hs_jonswap = 0.283 .* (tanh(0.53 .* g .* fetch ./ wind_speed.^2)).^0.44 .* wind_speed.^2 ./ g;
    Tp_jonswap = 0.84 .* (tanh(0.833 .* g .* fetch ./ wind_speed.^2)).^(-0.39) .* wind_speed ./ g;
    Tp_jonswap(isinf(Tp_jonswap)) = NaN;
    diff_H = abs(Hs_jonswap - Hs_model);
    diff_T = abs(Tp_jonswap - Tp_model);



    if make_plot

        disp('Significant Wave Heights Comparison:')
        disp(table(Hs_model(:), Hs_jonswap(:), 'VariableNames', {'Model','JONSWAP'}))
        fprintf('Max of %f m difference from JONSWAP\n',max(diff_H(:)));

        disp('Peak Wave Periods Comparison:')
        disp(table(Tp_model(:), Tp_jonswap(:), 'VariableNames', {'Model','JONSWAP'}))
        fprintf('Max of %f s difference from JONSWAP\n',max(diff_T(:)));

    end
end
