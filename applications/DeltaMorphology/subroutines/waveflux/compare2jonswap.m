function compare2jonswap(wind_speed, fetch, g, Hs_model, Tp_model)
%COMPARE2JONSWAP Compares model wave growth to JONSWAP growth curves



    % JONSWAP Growth Curves 
    Hs_jonswap = 0.283 .* (tanh(0.53 .* g .* fetch ./ wind_speed.^2)).^0.44 .* wind_speed.^2 ./ g;
    Tp_jonswap = 0.84 .* (tanh(0.833 .* g .* fetch ./ wind_speed.^2)).^(-0.39) .* wind_speed ./ g;
    Tp_jonswap(isinf(Tp_jonswap)) = NaN;

    % Differences
    diff_H = abs(Hs_model - Hs_jonswap);
    diff_T = abs(Tp_model - Tp_jonswap);

    
    figure('Name','Non-Dimensional Growth Curve Comparison');

    % Significant Wave Height Comparison

    clrs = winter(numel(wind_speed));
    
    subplot(2,1,1)
    for i = 1:numel(wind_speed)
        plot(fetch(i,:),Hs_model(i,:), '-', 'LineWidth', 2,'Color',clrs(i,:),'DisplayName',['u = ',num2str(wind_speed(i)),' Model']); 
        hold on;  
        plot(fetch(i,:),Hs_jonswap(i,:), '-o', 'LineWidth', 2,'Color',clrs(i,:),'DisplayName',['u = ',num2str(wind_speed(i)),' JONSWAP']);
    end
    grid on;
    xlabel('fetch');
    ylabel('H_s (m)');
    legend('show','Location','best');
    txt_title = sprintf('Hsig -- Max difference: %.3e m', max(diff_H(:)));
    title(txt_title);
    axis padded;

    subplot(2,1,2)
    for i = 1:numel(wind_speed)
        plot(fetch(i,:),Tp_model(i,:), '-', 'LineWidth', 2,'Color',clrs(i,:),'DisplayName',['u = ',num2str(wind_speed(i)),' Model']); 
        hold on;  
        plot(fetch(i,:),Tp_jonswap(i,:), '-o', 'LineWidth', 2,'Color',clrs(i,:),'DisplayName',['u = ',num2str(wind_speed(i)),' JONSWAP']);
    end
    grid on;
    xlabel('fetch');
    ylabel('H_s (m)');
    legend('show','Location','best');
    txt_title = sprintf('Tp -- Max difference: %.3e m', max(diff_T(:)));
    title(txt_title);
    axis padded;

end
