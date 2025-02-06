function sinuosity = calc_sinuosity(x,y,ws)

    % distance walking along the shore
    for i = 1:length(x)
        
        if i + 1 <= length(x)
            each_diff(i) = sqrt((y(i+1) - y(i))^2 + (x(i+1)-x(i))^2);
        else 
            each_diff(i) = sqrt((y(1) - y(i))^2 + (x(1)-x(i))^2);
        end
    end


   for i = 1:length(x)
    
        if i - ws >= 1 && i + ws <= length(x)
            j1 = i - ws;
            j2 = i + ws;
            full_length_alongshore(i) = sum(each_diff(j1:j2));
        elseif i - ws >= 1 && i + ws > length(x)
            j1 = i - ws;
            j2 = i + ws - length(x);
            full_length_alongshore(i) = sum(each_diff(j1:length(x))) + sum(each_diff(1:j2));      
        elseif i - ws < 1 && i + ws <= length(x)
            j1 = length(x) + (i - ws);
            j2 = i + ws;
            full_length_alongshore(i) = sum(each_diff(j1:length(x))) + sum(each_diff(1:j2));
        else
            disp('missed points')
            disp('i')
        end
        
        % P1 - > P2 is distance as crow flies 
        pt1 = [x(j1),y(j1)]; 
        pt2 = [x(j2),y(j2)];
        
        abs_diff(i) = sqrt((pt2(2) - pt1(2))^2 + (pt2(1) - pt1(1))^2); % absolute distance between two points on the shoreline          
           
   end

   % sinuosity is between 0 and 1
    %   s -> 1 means more sinous
    %   s -> 0 means less sinous
    sinuosity = abs_diff./full_length_alongshore;
    sinuosity = 1 - sinuosity;
    

end