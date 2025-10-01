clc
clear
close all

% estimate Q of Delta at Saraswati river mouth

g = 1.352;

Qs_river_labels = {
    'cold bedload min', 'cold bedload max';
    'warm bedload min', 'warm bedload max';
    'cold susload min', 'cold susload max';
    'warm susload min', 'warm susload max';
    };
% Birch+2023, appendix Table S5
Qs_river = [
    6.5e-4 0.36; % cold bedload
    8.9e-4 0.49; % warm bedload
    1.2 9.4      % cold suspended load
    0.43 3.3];   % warm suspended load


Qs_river = 0.00084978;

% Hayes+2010
A_slope = 2e-3; % section A from altimeter Hayes2010
L_slope = 1e-3; % section L from altimeter Hayes2010


H0 = 15;% 0.1; 
T = 1.5;% 6; 
L0 = (g*(T^2))/(2*pi);

rho = 588.15;
rho_s = 940;

D50 = 0.2;
nu = 8.084e-7;


alpha0 = deg2rad(0:90);
psi = deg2rad(0);
rel_angle = alpha0 - psi;
K = 0.77;%logspace(-2,2,10);
Qs_waves = calculate_volumetric_sed_flux_CERC(H0,T,alpha0,psi,K);


Qs_D = calculate_volumetric_sed_flux_Deigaard(H0,L0,rho_s,rho,g,D50,nu,slope);

figure;

% plot Sam's river input
for i = 1:numel(Qs_river)
    yline(Qs_river(i), '-', 'LineWidth', 3);
end
hold on

% colorcode the K value
nLines = 1;
colors = [linspace(0,1,nLines)', zeros(nLines,1), linspace(1,0,nLines)']; % blue → red

% Plot wave sed flux
for i = 1:nLines
    xvals = rad2deg(rel_angle);
    yvals = Qs_waves(i,:);
    
    plot(xvals, yvals, '--', 'LineWidth', 3, 'Color', colors(i,:));
    
    midIdx = round(numel(xvals)/2);
    x_text = xvals(midIdx);
    y_text = yvals(midIdx);
    
    if y_text > 1e-4
        text(x_text + 1, y_text* 10^(0.1), sprintf('%.2f', K(i)), ...
            'Color', colors(i,:), ...
            'FontSize', 10, ...
            'FontWeight', 'bold', ...
            'Clipping', 'on');  
    end
end


set(gca, 'YScale', 'log')
ylim([1e-5 1e5])
xlabel('relative approach angle (deg)')
ylabel('Qs (m^3/s)')
title('too small')

plot(rad2deg(alpha0),Qs_D,'LineWidth', 3, 'Color','b')
plot(rad2deg(alpha0),Qs_D + Qs_D.*0.5,'LineWidth',3,'Color','b')
plot(rad2deg(alpha0),Qs_D - Qs_D.*0.5,'LineWidth',3,'Color','b')
hold off;

nLines = numel(Qs_river);
colors = [linspace(0,1,nLines)', zeros(nLines,1), linspace(1,0,nLines)']; % blue → red

figure;

[Qs_river,ii] = sort(Qs_river(:));
Qs_river_labels = Qs_river_labels(:);
Qs_river_labels = Qs_river_labels(ii);
for i = 1:numel(Qs_river)
    R_delta(i,:) = Qs_river(i)./Qs_D;
    plot(rad2deg(alpha0),R_delta(i,:),'-','LineWidth',3,'Color',colors(i,:),'DisplayName',[Qs_river_labels{i},' : ',num2str(Qs_river(i))])
    hold on
end
legend('show')
set(gca, 'YScale', 'log')
yline(1,'--','HandleVisibility', 'off');
grid on
xlabel('relative approach angle (deg)')
ylabel('R')
title('too small')
% load('ontario_shoreline_pts.mat')
% lon_shore = mod(ontarioshorelinepts.x,360);
% lon_shore = [lon_shore ;lon_shore(1)];
% lat_shore = ontarioshorelinepts.y;
% lat_shore = [lat_shore ;lat_shore(1)];
% 
% 
% [x_shore, y_shore, distances,ref_pt] = titan_geodesic_distances(lat_shore,lon_shore);
% 
% % [A_Xmesh,A_Ymesh,A_zDep] = make_bathtub_lake(A_slope,[x_shore y_shore]);
% % [L_Xmesh,L_Ymesh,L_zDep] = make_bathtub_lake(L_slope,[x_shore y_shore]);
% load('A_slope.mat','A_Xmesh','A_Ymesh','A_zDep')
% load('L_slope.mat','L_Xmesh','L_Ymesh','L_zDep')
% 


function Qs = calculate_volumetric_sed_flux_Deigaard(H0,L0,rho_s,rho,g,D50,nu,alpha0,slope)

    w = calc_settling_velocity(rho_s,rho,g,D50,nu);
    R = rho_s/rho - 1;
    w_star = w/sqrt(g*D50);%(w.^3)./(R*g*nu);
    theta_0 = 0.1*((H0/D50)^2.3)*(sqrt(H0/L0))*exp(-6.1*w_star);
    
    alpha = alpha0./pi;
    theta = ((sin(2.*alpha0.*(1-(0.4.*(alpha).*(1-(alpha)))))).^(5/2)).*theta_0;
    s = rho_s/rho;
    Qs = theta*(H0*sqrt(slope)*sqrt(s-1)*g*D50);

    function settling_velocity = calc_settling_velocity(rho_s,rho,g,D50,nu)
        R = (rho_s/rho) - 1;
        D_star = (R*g.*(D50.^3))./(nu^2);
        log_D_star = log10(D_star);  
        log_W_star = -3.76715 + ...
                         1.92944 .* log_D_star - ...
                         0.09815 .* (log_D_star.^2) - ...
                         0.00575 .* (log_D_star.^3) + ...
                         0.00056 .* (log_D_star.^4);
        W_star = 10.^(log_W_star);
        settling_velocity = ((W_star.*R.*g.*nu)).^(1/3);

       
    end

end

function Qs = calculate_volumetric_sed_flux_CERC(H0,T,alpha0,psi,K)
% calculate CERC formula for volumetric sediment flux (m3/s)
% input:
%   H0      = deepwater wave height (m)
%   T       = deepwater wave period (s)
%   alpha0  = wave crest angle (rad)
%   psi     = shoreline angle (rad)
%   K       = empirical constant 
% output
%   Qs      = volumetric sediment flux (m3/s)

    
    rel_angle = alpha0 - psi;
    for i = 1:numel(K)
        Qs(i,:) = K(i)*(H0^(12/5))*(T^(1/5)).*(cos(rel_angle).^(6/5)).*sin(rel_angle);
    end

end

