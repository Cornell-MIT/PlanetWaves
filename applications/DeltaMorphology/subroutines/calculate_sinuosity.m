function S = calculate_sinuosity(x, y, windowsize)

x = x(:);
y = y(:);

n = length(x);
S = NaN(n,1);

halfWin = floor(windowsize/2);
windowsize = min(windowsize, n);

for i = 1:n

    idx = (i-halfWin):(i+halfWin);
    idx = mod(idx-1, n) + 1;

    xx = x(idx);
    yy = y(idx);

    % path length
    dpath = sum(sqrt(diff(xx).^2 + diff(yy).^2));

    % end-to-end distance
    dend = hypot(xx(end) - xx(1), yy(end) - yy(1));

    if dend > 0
        S(i) = dpath / dend;
    else
        S(i) = NaN;
    end
end

end