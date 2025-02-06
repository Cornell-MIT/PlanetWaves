function [shorelines_rotated] = rotate_basin(shoreline,hemisphere,plotting_on)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective:
%   This function will create the basin of interest and rotate so in sterographic projection but with (0,1) to be North pole
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
%   shoreline = shoreline coordinates for basin of interest  (in stereographic projection)
%               format: shorelines{1}(:,1) = x;  shorelines{1}(:,2) = y;
%   hemisphere = Northern (N) or Southern (S)
%   plotting_on = 1 (makes plots) else no plots printed
%
% outputs:
%   shorelines_rotated = coordinates from input shoreline rotated to right projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','all') % supress warnings from polyshape to printing to command line

disp('creating lake basin...')

x = shoreline{1}(:,1);
y = shoreline{1}(:,2);

polyin = polyshape({x},{y});
[xc,yc] = centroid(polyin);

% ROTATE MAP SO NORTH POLE IS IN +Y UNIT DIRECTION 

if hemisphere == 'N'
    theta_rot = -atan2(yc,xc) - pi/2; % in Northern Hemisphere (0,1) = North Pole
elseif hemisphere == 'S'
    theta_rot = -atan2(yc,xc) + pi/2; % in Southern Hemisphere, (0,-1) = South
else 
    disp('ERROR in make_basin -- hemisphere not specified')
end

v = [x y]';

center = repmat([xc; yc], 1, length(x));
R = [cos(theta_rot) -sin(theta_rot); sin(theta_rot) cos(theta_rot)];
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
x_rotated = vo(1,:);
y_rotated = vo(2,:);

dp0 = [0 0] - [xc yc];

s1 = [0 0]' - center;     % shift points in the plane so that the center of rotation is at the origin
so1 = R*s1;           % apply the rotation about the origin
vo1 = so1 + center;   % shift again so the origin goes back to the desired center of rotation
x0_rotated = vo1(1,end);
y0_rotated = vo1(2,end);

dp = [x0_rotated y0_rotated] - [xc yc];

if plotting_on
    f1 = figure(1);
    plot(x,y,'k')
    hold on
    plot(x(1:length(x)/5),y(1:length(x)/5),'g')
    plot(0,0,'xr') % the 0,0 point will rotate
    quiver(xc,yc,dp0(1),dp0(2),'AutoScale','off')
    hold off
    title('original basin')
    
    f2 = figure(2);
    plot(x_rotated,y_rotated','r')
    hold on
    plot(x_rotated(1:length(x)/5),y_rotated(1:length(x)/5),'g')
    plot(x0_rotated,y0_rotated,'xr') % rotated (0,0) point
    quiver(xc,yc,dp(1),dp(2),'AutoScale','off')
    title('rotated basin')
end

if hemisphere == 'N'
    disp('North pole in Northern Hemisphere should be at 90 degrees.')
elseif hemisphere == 'S'
    disp('North pole in Southern Hemisphere should be at -90 degrees.')
end
fprintf('North pole for rotated basin is at: %2.2f\n',rad2deg(atan2(dp(2),dp(1)))) % check that this is at 90 deg

shorelines_rotated{1}(:,1) = x_rotated;
shorelines_rotated{1}(:,2) = y_rotated;

disp('rotate_basin completed')
end