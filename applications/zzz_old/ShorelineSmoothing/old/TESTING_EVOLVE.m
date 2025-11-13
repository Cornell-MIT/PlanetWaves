clc
clear
close all


num_iterations = 100;
frac = [60 75 85 100];
u = [0 pi/2 pi 3*pi/2];
u = circshift(u,3);
alpha = 0.01;


num_points = 1000;
theta = linspace(0, 2*pi, num_points+1); 
theta(end) = [];
radius = 10;
x_smooth = radius * cos(theta);
x = x_smooth + 0.1.*randn(size(theta));
y_smooth = radius * sin(theta);
y = y_smooth + 0.1.*randn(size(theta)); 

x_original = x;
y_original = y;

[s_angle,relative_angle] = calc_shoreline_angle(x,y,0);
jet_wrap = vertcat(jet,flipud(jet));

figure; axis equal;
scatter(x, y, 50, rad2deg(relative_angle), 'filled'); 
colormap(jet_wrap)
colorbar

for i = 1:num_points
    % nieghbors
    prev_idx = mod(i-2, num_points) + 1;
    next_idx = mod(i, num_points) + 1;
    
    % avg position
    x_avg = (x(prev_idx) + x(next_idx)) / 2;
    y_avg = (y(prev_idx) + y(next_idx)) / 2;

end
H = ones(size(x));

figure;
xlim([-15 15])
ylim([-15 15])

samples = randperm(num_iterations);
for ii = 1:num_iterations
    
    my_sample = samples(ii);
    if my_sample <= frac(1)
        wind_direction = u(1);
        saved_wind(ii) = u(1);
    elseif my_sample > frac(1) && my_sample <= frac(2)
        wind_direction = u(2);
        saved_wind(ii) = u(2);
    elseif my_sample > frac(2) && my_sample <= frac(3)
        wind_direction = u(3);
        saved_wind(ii) = u(3);
    elseif my_sample > frac(3) && my_sample <= frac(4)
        wind_direction = u(4);
        saved_wind(ii) = u(4);
    else
        error('problem with sampler')
    end
    [x, y] = deform_shoreline(x, y, x_smooth,y_smooth, H, wind_direction, 1, alpha);
    plot(x_original,y_original,'--r')
    hold on
    plot(x,y,'-k')
    drawnow;
    hold off;
    pause(0.4)


end

figure;
plot([x_original x_original(1)],[y_original y_original(1)],'-r','LineWidth',2)
hold on
plot([x x(1)],[y y(1)],'-k')
legend('original','after diffusition')

sinuosity_before = calc_sinuosity(x_original,y_original,100);
sinuosity_after = calc_sinuosity(x,y,100);

figure;
plot(x_original, y_original, 'k-', 'LineWidth', 2);
scatter(x_original, y_original, 50, sinuosity_before, 'filled'); 
clim([1 10])
colorbar
title('sinuosity before diffusion')

figure;
plot(x, y, 'k-', 'LineWidth', 2);
scatter(x, y, 50, sinuosity_after, 'filled'); 
clim([1 10])
colorbar
title('sinuosity after diffusion')

figure;
histogram(sinuosity_before)
hold on
histogram(sinuosity_after)
legend('before','after')
xlabel('sinuosity')


function [shoreline_angle,relative_angle] = calc_shoreline_angle(x,y,wind_direction)
% Initialize arrays for theta and Qs, matching the order of shoreline points (x, y)
relative_angle = NaN(size(x));  % Relative angle of attack

% find theta
for i = 1:numel(x)
    if i < numel(x)
        shoreline_angle(i) = atan2(y(i+1) - y(i), x(i+1) - x(i));
    else % wrap around to the first point to complete the circle
        shoreline_angle(i) = atan2(y(1) - y(i), x(1) - x(i));
    end
    
    shoreline_angle(i) = mod(shoreline_angle(i),2*pi);
    wave_front_angle(i) =  wind_direction + pi / 2;
    relative_angle(i) = mod(wave_front_angle(i) - shoreline_angle(i),2*pi);


end

end

function [wave_energy,diffusivity] = calculate_wave_forcing(x,y,H,wind_direction)

    jet_wrap = vertcat(jet,flipud(jet));

    [s_angle,relative_angle] = calc_shoreline_angle(x,y,wind_direction);

    % figure; hold on; axis equal;
    % plot(x, y, 'k-', 'LineWidth', 2);
    % scatter(x, y, 50, rad2deg(s_angle), 'filled'); % Color by values
    % colormap(jet_wrap);
    % colorbar;
    % title('shoreline angle')
    % 
    % figure; hold on; axis equal;
    % plot(x, y, 'k-', 'LineWidth', 2);
    % quiver(0,0,5,0)
    % scatter(x, y, 50, rad2deg(relative_angle), 'filled'); % Color by values
    % colormap(jet_wrap);
    % colorbar;
    % title('angle relative to wave front angle')

    for i = 1:numel(x)
        wave_energy(i) = H(i)^2;
        wave_energy(i) = wave_energy(i)*cos(relative_angle(i));
        wave_energy(wave_energy<0) = 0;

        if relative_angle(i) >= 0 && relative_angle(i) <= pi/2
            diffusivity(i) = (H(i)^(12/5))*(cos(relative_angle(i))^(1/5))*(cos(relative_angle(i))^2 - (6/5)*sin(relative_angle(i))^2);
        else
            diffusivity(i) = 0;
        end
        % if ~isreal(diffusivity(i))
        %     diffusivity(i) = 0;
        % end
    end

end

function [x_out, y_out] = deform_shoreline(x, y, x_avg,y_avg, H, wind_direction, num_iterations, alpha)

    num_points = length(x);

    for iter = 1:num_iterations
        x_new = x;
        y_new = y;

       % calculate diff
        [~,z] = calculate_wave_forcing(x,y,H,wind_direction);

        for i = 1:num_points
            % move relative to avg depending on diffusivity
            if z(i) > 0
                x_new(i) = x(i) + alpha*z(i)*(x_avg(i) - x(i));
                y_new(i) = y(i) + alpha*z(i)*(y_avg(i) - y(i));
            elseif z(i) < 0
                x_new(i) = x(i) - alpha*abs(z(i))*(x_avg(i) - x(i));
                y_new(i) = y(i) - alpha*abs(z(i))*(y_avg(i) - x(i));
            end
         
        end

        % Update positions
        x = x_new;
        y = y_new;

    end

    % Output final coordinates
    x_out = x;
    y_out = y;

end