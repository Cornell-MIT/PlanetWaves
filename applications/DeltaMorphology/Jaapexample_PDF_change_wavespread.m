clc;
clear;
close all

% Nienhuis 2015, supp fig A,C
H0 = 1;
T0 = 1;
main_wind_dir = 0;
kappa = 0.1:0.5:10;
kappa = flip(kappa); 
theta_rad = deg2rad(wrap_180(rad2deg(0:0.01*pi:2*pi)));  
theta_deg = rad2deg(theta_rad);

phi0 = deg2rad(main_wind_dir);

LST = nan(numel(kappa), numel(theta_rad));
Qs = nan(numel(kappa), numel(theta_rad));
E_pdf = nan(numel(kappa), numel(theta_rad));
integral_Qs = nan(size(kappa));

gray_colors = repmat(linspace(0, 0.8, numel(kappa))', 1, 3); % Each row is [R G B]

fig = figure('Position', [100 100 900 600]);
tiledlayout(2,1)

ax1 = nexttile;
hold on
for m = 1:numel(kappa)
    E_pdf(m,:) = vonMises(phi0, kappa(m), theta_rad) ;

    plot(rad2deg(phi0 - theta_rad), E_pdf(m,:), 'Color', gray_colors(m,:), 'LineWidth', 1.5)
end
ylabel('E pdf')
xlabel('phi0 - theta')
grid on
xlim([-200 200])

ax2 = nexttile;
hold on
for m = 1:numel(kappa)
    relative_angle = deg2rad(wrap_180(rad2deg(phi0 - theta_rad)));

    LST(m,:) = CERC(H0, T0, relative_angle);
    

    Qs(m,:) = E_pdf(m,:) .* LST(m,:);
    integral_Qs(m) = round(sum(abs(Qs(m,:)) * (theta_rad(2) - theta_rad(1)),'omitmissing'),5);

    [minQ,iminQ] = min(Qs(m,:));
    [maxQ,imaxQ] = max(Qs(m,:));

    plot(theta_deg, Qs(m,:), 'Color', gray_colors(m,:), 'LineWidth', 1.5)
    plot(theta_deg(iminQ),minQ,'ob')
    plot(theta_deg(imaxQ),maxQ,'or')
end
xlabel('\theta (deg)')
ylabel('Qs')
grid on
xlim([-200 200])

figure;
plot(kappa,integral_Qs,'ok','LineWidth',3)
hold on;
xlabel('kappa')
ylabel('integral of Qs over all shoreline angles theta')

[minWAVE,iminWAVE] = min(integral_Qs);
[maxWAVE,imaxWAVE] = max(integral_Qs);

plot(kappa(iminWAVE),minWAVE,'ob','MarkerSize',10,'MarkerFaceColor','b')
plot(kappa(imaxWAVE),maxWAVE,'or','MarkerSize',10,'MarkerFaceColor','r')

