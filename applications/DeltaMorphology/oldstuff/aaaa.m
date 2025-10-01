clc;
clear;
close all

% Parameters
H0 = 1;
T0 = 1;
main_wind_dir = 0;
kappa = 0.1:0.5:10;
kappa = flip(kappa); % High kappa first = black lines
theta_deg = linspace(-180, 179, 359);  
theta_rad = deg2rad(theta_deg);

phi0 = deg2rad(main_wind_dir);

% Preallocate
LST = nan(numel(kappa), numel(theta_rad));
Qs = nan(numel(kappa), numel(theta_rad));
E_pdf = nan(numel(kappa), numel(theta_rad));
integral_Qs = nan(size(kappa));

% Grayscale colormap: black to light gray
gray_colors = repmat(linspace(0, 0.8, numel(kappa))', 1, 3); % Each row is [R G B]

% Create figure and subplots
fig = figure('Position', [100 100 900 600]);
tiledlayout(2,1)

% === Top subplot: von Mises PDF ===
ax1 = nexttile;
hold on
for m = 1:numel(kappa)
    % von Mises PDF
    E_pdf(m,:) = vonMises(phi0, kappa(m), theta_rad) ;

    % Plot E_pdf with corresponding grayscale
    plot(phi0 - theta_deg, E_pdf(m,:), 'Color', gray_colors(m,:), 'LineWidth', 1.5)
end
ylabel('E pdf')
xlabel('phi0 - theta')
grid on
xlim([-200 200])

% === Bottom subplot: Qs ===
ax2 = nexttile;
hold on
for m = 1:numel(kappa)
    % Relative angles
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
    Qs(m,:) = E_pdf(m,:) .* LST(m,:);
    integral_Qs(m) = sum(Qs(m,:) * (theta_rad(2) - theta_rad(1)));

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
