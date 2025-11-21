function shoreline_angle = calc_regional_shoreline_angle(x, y,window_size)
% calculates the regional shoreline angle for the shoreline points (x,y)
% over a window size window_size. Possible shoreline angles between 
% -180 to 180 deg (atan2), defined CCW.
% E.g.,
% shoreline_angle = 0 (shoreline goes from left to the right)
% shoreline_angle = pi/2 (shoreline goes from bottom to top)
% shoreline_angle = pi (shoreline goes from right to left)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%   x             = x-coordinates of shoreline [1 x numel(points)] 
%   y             = y-coordinates of shoreline [1 x numel(points)]
%   window_size   = number of points computing the angle for chord distance over [1 x 1]
% OUTPUTS:
% shoreline_angle = regional shoreline angle tangential to shoreline [rad] [1 x numel(points)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    make_plot = 1; % turn on to make a plot of the vector defined by the 
                   % regional shoreline angle (tangential to shore)
    
    % force CCW so normal points in to lake and not out from lake
    [x,y] = make_CCW(x,y);
    
    shoreline_angle = NaN(1,numel(x));
    
    for i = 1:length(x)-1 % ignore the repeated element at the end
        
        % circular forward index 
        idx_f = i + window_size;
        if idx_f > length(x)
            idx_f = idx_f - length(x);
        end
        
        % circular backward index 
        idx_b = i - window_size;
        if idx_b < 1
            idx_b = idx_b + length(x);
        end
        
        % chord direction across the window size
        dx = x(idx_f) - x(idx_b);
        dy = y(idx_f) - y(idx_b);
        
        shoreline_angle(i) = atan2(dy, dx);
    end
   
    shoreline_angle(end) = shoreline_angle(1); % because closed shape
    
    if make_plot
           figure('Name','Regional Shoreline Tangential Vector'); 

           plot(x, y, 'k-', 'LineWidth', 1.5);
           hold on
           fill(x,y,[0 0.7 0.7],'FaceAlpha',0.2)
           scatter(x,y,50,rad2deg((shoreline_angle)),'filled')
           cb = colorbar();
           ylabel(cb,'Shoreline Angle [Deg]','FontSize',16, 'Rotation',270)
           scale_arrow = 1;
           u = cos(shoreline_angle);
           v = sin(shoreline_angle);
           quiver(x, y, scale_arrow.*u, scale_arrow.*v, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
           axis equal padded
    end
   
end