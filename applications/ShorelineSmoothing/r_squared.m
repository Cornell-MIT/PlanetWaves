function [R2,slope] = r_squared(x, y,xval,yval)
    x = x(:); % Ensure x is a column vector
    y = y(:); % Ensure y is a column vector

    % Remove NaN values in x or y
    validIdx = ~isnan(x) & ~isnan(y);
    x = x(validIdx);
    y = y(validIdx);

    p = polyfit(x, y, 1); % Fit a linear model
    slope = p(1);
    y_fit = polyval(p, x); % Compute fitted values
    
    SS_tot = sum((y - mean(y)).^2); % Total sum of squares
    SS_res = sum((y - y_fit).^2); % Residual sum of squares
    
    R2 = 1 - (SS_res / SS_tot); % Compute R^2 value
    
    % Plot the data and the linear fit
    figure;
    plot(x, y, 'bo', 'MarkerFaceColor', 'b'); % Original data points
    hold on;
    plot(x, y_fit, 'r-', 'LineWidth', 2); % Fitted line
    xlabel(xval);
    ylabel(yval);
    title('Linear Fit');
    legend('Data', 'Fit');
    grid on;
    hold off;
    title(sprintf('R^2 = %f',R2))
end