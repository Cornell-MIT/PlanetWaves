clc;
clear;
close all;

make_gif = 1;
% Nienhuis 2015, supp fig B,D

% shoreline orientation
my_shoreline = deg2rad([90:-10:-90]);
% DOMINANT WIND DIRECTION(S)
phi0 = deg2rad(0);  
% ALL POSSIBLE RELATIVE ORIENTATIONS
rel_angle = deg2rad(-89:89);  

% CONSTANTS FOR WAVE FIELD
H0 = 1;
T0 = 1;
kappa = 10;  

original_shoreline = deg2rad(360);
figure;
quiver(0,0,cos(phi0),sin(phi0),'k');
hold on
arrow_length = 2;
x0 = 0;
y0 = -1;
dx = arrow_length * cos(original_shoreline);
dy = arrow_length * sin(original_shoreline);
quiver(x0, y0, dx, dy, 'r')

% Compute perpendicular direction (normal vector)
% Rotate 90 degrees clockwise for 'land' (right side)
% Rotate 90 degrees counter-clockwise for 'sea' (left side)
normal_length = 1; % distance from arrow for label
nx = normal_length * cos(original_shoreline - pi/2); % right side
ny = normal_length * sin(original_shoreline - pi/2);
lx = normal_length * cos(original_shoreline + pi/2); % left side
ly = normal_length * sin(original_shoreline + pi/2);

% Position labels
text(x0 + dx/2 + nx, y0 + dy/2 + ny, 'land')
text(x0 + dx/2 + lx, y0 + dy/2 + ly, 'sea')

% Axis setup
axis equal
axis([-3 3 -3 3])

% Define grey color for past plots
grey_color = [0.5, 0.5, 0.5]; % Medium grey color for past plots

% Initialize cell arrays for plot handles
wave_energy_handles = {};  
sediment_transport_handles = {};  % Use cell array to store the handles

% Compute Wave Energy PDF and Sediment Transport Function
for m = 1%:numel(my_shoreline)

    % Calculate relative wind direction
    wind = rel_angle - my_shoreline(m);  % Relative wind direction to shoreline

    % Filter valid wind directions within range [-pi/2, pi/2]
    valid_wind = (rad2deg(wind) > 90);

    % Von Mises PDF (you may need to define the vonMises function separately)
    E_pdf(m,:) = vonMises(phi0 - my_shoreline(m), kappa, wind);  % Energy PDF (Von Mises)
    E_pdf(m,:) = delta_function(E_pdf(m,:));

    % Compute CERC (sediment transport function, assuming it's a known function)
    LST = CERC(H0, T0, rel_angle);  % LST(wind-shoreline)

    % Sediment transport (Qs)
    Qs = E_pdf(m,:) .* LST;  % Combined sediment transport

    theta =  rel_angle + my_shoreline(m);
    valid_theta = (theta >= -pi/2) & (theta <= pi/2);
    Qs(~valid_theta) = NaN;

    % Find min and max Qs
    [minQ(m), iminQ(m)] = min(Qs);
    [maxQ(m), imaxQ(m)] = max(Qs);

    sum_Qs(m) = maxQ(m) - minQ(m);  % Net transport max

    if make_gif
        % Plot Wave Energy
        subplot(2,1,1);
        title(sprintf('For a shoreline $\theta_0$ with original %.0f deg orientation',rad2deg(my_shoreline(m))))
    
        if m == 1
            % Plot the first curve in black
            h1 = plot(rad2deg(wind), E_pdf, 'k', 'LineWidth', 1.5);  
        else
            % Update previous curve color to grey and plot new curve in black
            for k = 1:length(wave_energy_handles)
                set(wave_energy_handles{k}, 'Color', grey_color);  % Set previous curves to grey
            end
            h1 = plot(rad2deg(wind), E_pdf, 'k', 'LineWidth', 1.5);  % Plot current curve in black
        end
        wave_energy_handles{end+1} = h1;  % Store the current handle for future updates
    
        hold on;
        xlim([-200 200])
        ylim([0 2])
        xlabel('wave angle \phi_0 (deg)')
        ylabel('wave energy PDF')
    
        % Plot Sediment Transport
        subplot(2,1,2); 
        hold on;
        if m == 1
            % Plot the first curve in black
            h2 = plot(rad2deg(theta), Qs, 'k', 'LineWidth', 1.5);  
        else
            % Update previous curve color to grey and plot new curve in black
            for k = 1:length(sediment_transport_handles)
                set(sediment_transport_handles{k}, 'Color', grey_color);  % Set previous curves to grey
            end
            h2 = plot(rad2deg(theta), Qs, 'k', 'LineWidth', 1.5);  % Plot current curve in black
        end
        % Plot min and max Qs values
        plot(rad2deg(theta(iminQ(m))), minQ(m), 'ob')
        plot(rad2deg(theta(imaxQ(m))), maxQ(m), 'or')
        xlim([-100 100])
        ylim([-0.3 0.3])
        xlabel('delta shoreline orientation \theta (deg)')
    
        sediment_transport_handles{end+1} = h2;  % Store the current handle for future updates
    
        if m == 1
            gif('movingshoreling.gif','DelayTime',0.5,'overwrite',true)
        else
            gif;
        end
    end
end

figure;
plot(rad2deg(my_shoreline), round(sum_Qs, 5), '-ok', 'MarkerFaceColor', 'k')
xlabel('shoreline orientation \theta_0 (deg)')
ylabel('max Qs left + max Qs right')


% figure;
% 
% 
% phi0 = 0;
% theta = deg2rad(-90:90);
% theta = flip(theta);
% 
% for i = 1:numel(theta)
% 
%     subplot(2,1,1);
% 
%     quiver(0,0,-4,0,'LineWidth',1.5,'MaxHeadSize',10)
%     hold on
% 
%     sx = -15*cos(theta(i));
%     sy = 15*sin(theta(i));
%     quiver(0,0,sx,sy,'g','LineWidth',1.5)
%     quiver(0,0,-sx,-sy,'g','ShowArrowHead','off','LineWidth',1.5)
% 
% 
% 
%     v_len = sqrt(sx^2 + sy^2);
%     ux = -sy / v_len; 
%     uy =  sx / v_len;  
% 
%     % Length of small arrow
%     arrow_len = 3;
% 
%     % Start point: tip of main arrow
%     start_x = sx - 0.1*sx;
%     start_y = sy - 0.1*sy;
% 
%     quiver(start_x, start_y, arrow_len*ux, arrow_len*uy,'r', 'LineWidth', 1.2, 'MaxHeadSize', 2);  
%     axis([-20 20 -20 20])
%     hold off;
% 
% 
%     subplot(2,1,2); 
%     plot(rad2deg(phi0 - theta),zeros(size(theta)),'-k')
%     hold on
%     plot(rad2deg(phi0 - theta(i)),0,'or','MarkerFaceColor','r')
%     hold off
%     axis([-100 100 -1 1])
%     xlabel('\phi_0 - \theta')
%     drawnow
% 
%     if i == 1
%         gif('shorelineorientation.gif')
%     else
%         gif;
%     end
% 
% 
% end



function simple_array = delta_function(original_array)

    [~,i] = max(original_array,[],"omitmissing");
    simple_array = zeros(size(original_array));
    simple_array(i) = 1;
    

end