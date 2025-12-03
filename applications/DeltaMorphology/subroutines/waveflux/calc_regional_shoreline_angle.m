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
                
        A = i - window_size; % start of window
        B = i + window_size; % end of window
        if A < 1
            A = length(x) + A;
        end
        if B > length(x)
            B = B - length(x);
        end
        
        
        % chord direction across the window size
        dx = x(B) - x(A);
        dy = y(B) - y(A);

        
        
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
           L = sqrt(dx.^2 + dy.^2);
           u = cos(shoreline_angle);
           v = sin(shoreline_angle);
           quiver(x, y, L .* u, L .* v, 0, 'r', 'LineWidth', 2, 'ShowArrowHead', 'off');
           axis equal padded
    end
   
end