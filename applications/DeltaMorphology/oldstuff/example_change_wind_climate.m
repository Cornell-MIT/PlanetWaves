
clc;
clear;
close all

% Parameters
H0 = 1;
T0 = 1;
kappa = 20;  
main_wind_dir_deg = 10;  

theta_deg = linspace(-180, 179, 359);  
theta_rad = deg2rad(theta_deg);
num_dirs = numel(main_wind_dir_deg);

LST = nan(num_dirs, numel(theta_rad));
Qs = nan(num_dirs, numel(theta_rad));
E_pdf = nan(num_dirs, numel(theta_rad));
integral_Qs = nan(size(main_wind_dir_deg));

gray_colors = repmat(linspace(0, 0.8, num_dirs)', 1, 3);  

fig = figure('Position', [100 100 900 600]);
tiledlayout(2,1)

ax1 = nexttile;
hold on
for m = 1:num_dirs
    phi0 = deg2rad(wrap_180(main_wind_dir_deg(m)));

    % von Mises PDF
    E_pdf(m,:) = vonMises(phi0, kappa, theta_rad);


    % Plot E_pdf
    plot(phi0 - theta_deg, E_pdf(m,:), 'Color', gray_colors(m,:), 'LineWidth', 1.5)
    xlabel('phi0 - theta')
    fprintf('E_pdf max-min at m=%d: %.5f\n', m, max(E_pdf(m,:)) - min(E_pdf(m,:)));

end
ylabel('E pdf')
grid on
xlim([-180 180])

ax2 = nexttile;
hold on
for m = 1:num_dirs
    phi0 = deg2rad(wrap_180(main_wind_dir_deg(m)));
    relative_angle = phi0 - theta_rad;

    % Compute LST
    for i = 1:numel(theta_rad)
        if relative_angle(i) > -pi/2 && relative_angle(i) < pi/2
            LST(m,i) = CERC(H0, T0, relative_angle(i));
        else
            LST(m,i) = 0;
        end
    end


  % Compute Qs
Qs_calc = E_pdf(m,:) .* LST(m,:);

% Define valid angle range
valid_idx = theta_rad > -pi/2 & theta_rad < pi/2;
Qs_calc(~valid_idx) = 0;  % Use NaN so they are ignored in calculations

Qs(m,:) = Qs_calc;

% Integrate Qs
Qs_summed(m) = trapz(theta_rad, Qs(m,:));

% Plot Qs
plot(theta_deg, Qs(m,:), 'Color', gray_colors(m,:), 'LineWidth', 1.5); hold on;

    



end
xlabel('\theta (deg)')
ylabel('Qs')
grid on


figure;
plot(main_wind_dir_deg,Qs_summed)
xlabel('\phi0 [deg]')
ylabel('sum Qs over all shoreline directions')