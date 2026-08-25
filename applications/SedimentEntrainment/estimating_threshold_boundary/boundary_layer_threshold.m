clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the threshold for oscillatory flow entrainment given spread in
% data as seen in FIGURE 4 (Komar and Miller 1973). This data was digitized
% and then fit with a line below and above the cutoff for laminar or 
% turbulent particle Reynolds numbers. The digitized data is stored in
% lam_turb_data.csv. The fits for the thresholds will be reported to the
% console.
% PLOTS:
%   (1) ratio of grain size and orbital diameter size to  particle Re
% Author: Una Schneck
% Last Modified: 8/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load in data from Komar and Miller 1973
data = readtable('lam_turb_data.csv');
% X = sqrt(d0/D)
% Y = rho * um^2 / ((rho_s - rho)gD)

% seperate laminar and turbulent data 
X_laminar = data.x(data.y>12);
Y_laminar = data.y(data.y>12);
X_turb = data.x(data.y<=12);
Y_turb = data.y(data.y<=12);

% plot the data
figure('Name','Komar & Miller 1973; Figure 4');
plot(data.x,data.y,'ok')
hold on;
plot(X_laminar,Y_laminar,'ok','MarkerFaceColor','r')
plot(X_turb,Y_turb,'ok','MarkerFaceColor','b')

% laminar_threshold_manohar_reported = 0.392.*X_laminar;
% plot(X_laminar,laminar_threshold_manohar_reported,'-r','LineWidth',3)

% fit the data for turbulent and laminar seperatly
coef_turb =  polyfit(X_turb,Y_turb,1);
turbulent_threshold_manohar = polyval(coef_turb,X_turb);

coef_lam = polyfit(X_laminar,Y_laminar,1);
laminar_threshold_manohar = polyval(coef_lam,X_laminar);

% plot the linear fits
plot(X_turb,turbulent_threshold_manohar,'-b','LineWidth',3)
plot(X_laminar,laminar_threshold_manohar,'-r','LineWidth',3)
grid on;


%plot(X_laminar,0.3*X_laminar,'-k','LineWidth',3) % this is the final reported threshold in the paper
% coef2 = polyfit(X_laminar,Y_laminar,1);
% lam2 = polyval(coef2,X_laminar);
% coef3 = polyfit(data.x,data.y,1);
% all_fit = polyval(coef3,data.x);
%plot(data.x,lam2,'--k')
%plot(data.x,all_fit,':k')

title('Komar & Miller 1973; Figure 4 (Manohar Data)')

yline(12,':k','aprox cutoff for transition from smooth to rough boundary')

xlabel('$(d_0/D_{50})^{0.5}$','Interpreter','latex')
ylabel('$\rho u_m^2 / ((\rho_s - \rho)gD_{50})$','Interpreter','latex')

disp('===========')
fprintf('Laminar Threshold : %f sqrt(d0/D50)\n',coef_lam(1))
fprintf('Turbulent Threshold : %f sqrt(d0/D50)\n',coef_turb(1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAME DATA IN FIGURE 5 AS FIGURE 4 FOR TURBULENT BOUNDARY LAYER POINTS BUT USING A DIFFERENT
% Y AXIS THAT CANT BE CONVERTED WITHOUT T AND um

% FIGURE 5 (Komar and Miller 1973)
% find approx fit to Manohar 1955 data for turbulent boundary layer (Fig 5, Komar & Miller 1973)
% X = d0/D
% Y = (rho*um)/(rho_s - rho)gT
%trend = readtable('manohar_1955_digitized_fig5.csv');
% X = trend.x;
%
% figure;
% plot(X,trend.y,'or','MarkerFaceColor','r','DisplayName','Komar & Miller Curve')
% 
% hold on
% 
% 
% fitnum = 3;
% coef = polyfit(log10(trend.x),log10(trend.y),fitnum);
% y_fit = polyval(coef,log10(trend.x));
% y_fit = 10.^y_fit;
% plot(X,y_fit,'-k','DisplayName',['Polyfit: n = ',num2str(fitnum)])
% 
% 
% laminar_threshold = 0.3.*sqrt(trend.x);
% plot(X,laminar_threshold,'--k','DisplayName','laminar threshold')
% 
% title(sprintf('Y = 10^{%f * (X^3) + %f * (X^2) + %f * (X) + %f}',coef(1),coef(2),coef(3),coef(4)))
% xlabel('X = $d_0/D$','Interpreter','latex')
% ylabel('Y  = $(\rho*u_m)/(\rho_s - \rho)gT$','Interpreter','latex')
% grid on;
% 
% legend('show','Location','best')
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% % xlim([100 3000])
% % ylim([0.003 0.03])