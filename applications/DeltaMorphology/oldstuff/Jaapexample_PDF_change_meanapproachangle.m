clc;
clear;
close all

% Nienhuis 2015, supp fig B,D

% CONSTANTS FOR WAVE FIELD
H0 = 1;
T0 = 1;
kappa = 10; 

% DOMINANT WIND DIRECTION(S)
main_wind_dir_rad = deg2rad(wrap_180(rad2deg(0:0.01*pi:2*pi)));  
main_wind_dir_deg = rad2deg(main_wind_dir_rad);

% ALL SHORELINE ORIENTATIONS
theta_rad = deg2rad(wrap_180(rad2deg(0:0.01*pi:2*pi)));  
theta_deg = rad2deg(theta_rad);


LST = nan(numel(main_wind_dir_rad), numel(theta_rad));
Qs = nan(numel(main_wind_dir_rad), numel(theta_rad));
E_pdf = nan(numel(main_wind_dir_rad), numel(theta_rad));
integral_Qs = nan(size(main_wind_dir_rad));
gray_colors = repmat(linspace(0, 0.8, numel(main_wind_dir_rad))', 1, 3); 

fig = figure('Position', [100 100 900 600]);
tiledlayout(2,1)

% COMPUTE WAVE ENERGY PDF
ax1 = nexttile;
hold on
for m = 1:numel(main_wind_dir_rad)
    phi0 = main_wind_dir_rad(m); 
    E_pdf(m,:) = vonMises(phi0, kappa, theta_rad);

    plot(phi0 - theta_deg, E_pdf(m,:), 'Color', gray_colors(m,:), 'LineWidth', 1.5)
end
ylabel('E pdf')
xlabel('phi0 - theta')
grid on
xlim([-200 200])

ax2 = nexttile;
hold on

% Calculate LST and Qs
for m = 1:numel(main_wind_dir_rad)
    phi0 = main_wind_dir_rad(m); 
    relative_angle = phi0 - theta_rad;

    % Compute LST
    LST(m,:) = CERC(H0, T0, relative_angle);

% INTEGRATE QS OVER ALL SHORELINE ORIENTATIONS
    Qs(m,:) = E_pdf(m,:) .* LST(m,:);
    integral_Qs(m) = round(sum(abs(Qs(m,:)) * (theta_rad(2) - theta_rad(1)),'omitmissing'),5);

    [minQ,iminQ] = min(Qs(m,:));
    [maxQ,imaxQ] = max(Qs(m,:));

    % Plot Qs with same grayscale
    plot(theta_deg, Qs(m,:), 'Color', gray_colors(m,:), 'LineWidth', 1.5)
    plot(theta_deg(iminQ),minQ,'ob')
    plot(theta_deg(imaxQ),maxQ,'or')
end
xlabel('\theta (deg)')
ylabel('Qs')
grid on
xlim([-200 200])


figure;
plot(main_wind_dir_deg,integral_Qs,'ok','LineWidth',3)
hold on;
xlabel('phi0')
ylabel('integral of Qs over all shoreline angles theta')

[minWAVE,iminWAVE] = min(integral_Qs);
[maxWAVE,imaxWAVE] = max(integral_Qs);

plot(main_wind_dir_deg(iminWAVE),minWAVE,'ob','MarkerSize',10,'MarkerFaceColor','b')
plot(main_wind_dir_deg(imaxWAVE),maxWAVE,'or','MarkerSize',10,'MarkerFaceColor','r')

