function [xnew,ynew] = smooth_path(x,y,n)

numPoints = floor(length(x) / n);

xnew = zeros(1, numPoints);
ynew = zeros(1, numPoints);

for i = 1:numPoints
    startIdx = (i - 1) * n + 1;
    endIdx = startIdx + n - 1;
    xnew(i) = mean(x(startIdx:endIdx),"omitmissing");
    ynew(i) = mean(y(startIdx:endIdx),"omitmissing");
end

if xnew(1) ~= xnew(end) || ynew(1) ~= ynew(end)
    xnew = [xnew xnew(1)];
    ynew = [ynew ynew(1)];
end
end

