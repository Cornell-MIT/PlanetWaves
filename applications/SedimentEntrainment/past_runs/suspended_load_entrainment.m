clc
clear
close all


% g = 9.81;
% rho_s = 2000;
% rho = 998.21;
% R = (rho_s/rho) - 1;
% H0 = 0:0.1:6;
% T = linspace(0,5,numel(H0));
% T(T==0) = NaN;

%load('C:\Users\Owner\OneDrive\Documents\00_Main\Work\Github_Repos\PlanetWaves\applications\SedimentEntrainment\past_runs\Waves03_43.mat','Planet','H0','T')
%g = Planet.gravity;

d50 = [6.35e-5 0.1]; % grain size


pp = 1;

planet_name = {'Earth','Titan'};
planet_g = [9.81 1.352];
planet_i = [1 4];
planet_rho = [1000 600];
planet_rho_s = [2000 1000];
planet_nu = [1.2e-6 8e-7];

load('AllPlanets_together_longerT.mat','htgrid','T_p','u_of_planet')
g = planet_g(pp);
m = planet_i(pp);
rho = planet_rho(pp);
rho_s = planet_rho_s(pp);
nu = planet_nu(pp);

winds = u_of_planet{pp};
max_wind = find(winds==5);


for i =  1:max_wind
    
    if ~isempty(htgrid{m,i})
        H0(i) = htgrid{m,i}{end}(5,5);
        T(i) = T_p{m,i}(5,5);
    else
        H0(i) = 0;
        T(i) = NaN;
    end

end

for comp = 1
    for i = 1:numel(H0)
        if H0(i) > 0 
        [~, Hshoal, ~, ubshoal,Ashoal,~] = shoaling_to_break_slope(H0(i), T(i), 80, 1e-3,g);
            [Hmax,iHmax] = max(Hshoal);
            a(comp,i) = Hmax;
            ub(comp,i) = ubshoal(iHmax);
            Ab(comp,i) = Ashoal(iHmax);
        else
            a(comp,i) = NaN;
            ub(comp,i) = NaN;
            Ab(comp,i) = NaN;
        end
    end
end


figure;
plot(winds(1:max_wind),Ab,'LineWidth',3)
xlabel('wind speed')
ylabel('max orbital size')
grid on
title(planet_name{pp})


% delta_ks = logspace(-2,2,100); % > 0.25 roughness can be felt
% fw = solve_for_fw(delta_ks);

delta_ks_sand = solve_for_delta(Ab,d50(1))/d50(1);                         % boundary layer thickness for sand / Nikuradse roughness for sand (sand size)
fw_sand = solve_for_fw(delta_ks_sand);                                     % wave friction factor for sand
delta_ks_grav = solve_for_delta(Ab,d50(2))./d50(2);                        % boundary layer thickness for gravel / Nikuradse roughness for gravel (gravel size)
fw_grav = solve_for_fw(delta_ks_grav);                                     % wave friction factor for gravel

figure
hold on
plot(delta_ks_sand,fw_sand,'-r','LineWidth',5,'DisplayName',['k_s = ', num2str(d50(1)), ' m'])
plot(delta_ks_grav,fw_grav,'-b','LineWidth',5,'DisplayName',['k_s = ', num2str(d50(2)), ' m'])
xlabel('ratio of boundary layer thickness and Nikuradse roughness $\delta/k_s$','Interpreter','latex')
ylabel('wave friction factor f$_w$','Interpreter','latex')
grid on
legend('show')
title(planet_name{pp})

d50_range = logspace(-5,-1,100);
Re_p = calc_particle_Reynolds(rho_s,rho,g,d50_range,nu);
threshold_entrainment1 = calculate_Parker_shields_curve(Re_p);
[suspension,s2,s3,s4] = calc_suspension_threshold(rho_s,rho,g,d50_range,nu);

mask_valid = suspension > threshold_entrainment1;

figure;
plot(Re_p,threshold_entrainment1,'-k')
hold on
plot(Re_p(mask_valid),suspension(mask_valid),'-k')
set(gca,'XScale','log')
set(gca,'YScale','log')



min_depth = 10^-4;
lake_slope = 2e-3;
d = 80:-lake_slope:min_depth;

for z = numel(d)

    depth = 1;
    for i = 1:numel(ub)
    
        Re_p_sand(i) = calc_particle_Reynolds(rho_s,rho,g,d50(1),nu);
        Re_p_grav(i) = calc_particle_Reynolds(rho_s,rho,g,d50(2),nu);
        theta_sand(i) = calc_shields_number(rho_s,rho,depth,ub(i),d50(1),g,min(fw_sand));
        theta_grav(i) = calc_shields_number(rho_s,rho,depth,ub(i),d50(2),g,min(fw_grav));
    
    end
    
   
    plot(Re_p_sand,theta_sand,'or','MarkerFaceColor','r')
    plot(Re_p_grav,theta_grav,'ob','MarkerFaceColor','b')
    drawnow
    title(['Depth = ',num2str(depth),' m (Distance ~',num2str(depth/lake_slope),' m offshore)'])
end
grid on
%axis equal
axis([10^-3 10^6 10^-4 10^6])
xlabel('Re$_p$','Interpreter','latex','FontSize',15)
ylabel('$\theta$','Interpreter','latex','FontSize',15)

H_plot_Titan = 0.01:0.1:3;
H_plot_Earth = 0.01:0.01:0.4;

figure;
plot(H_plot_Titan,solve_for_delta(H_plot_Titan,d50(1))./d50(1),'-r','LineWidth',3,'DisplayName','Titan, Sand')
hold on
plot(H_plot_Titan,solve_for_delta(H_plot_Titan,d50(2))./d50(2),'-b','LineWidth',3,'DisplayName','Titan, Gravel')
plot(H_plot_Earth,solve_for_delta(H_plot_Earth,d50(1))./d50(1),':r','LineWidth',4,'DisplayName','Earth, Sand')
plot(H_plot_Earth,solve_for_delta(H_plot_Earth,d50(2))./d50(2),':b','LineWidth',4,'DisplayName','Earth, Gravel')
xlabel('Wave Orbital Size')
ylabel('$\frac{\delta}{D_{50}}$','Interpreter','latex','FontSize',15)
legend('show')
grid on

function delta = solve_for_delta(H, K)
% Solve 30*(delta/k)*log10(30*delta/K) = 1.2*H/K for delta
% H > 0, H/K >> 1

    converged = 0;
    H = H(:);

    % Initial guess is asymptotic formula
    Z = (1.2*H*log(10))/K;
    delta = (0.04*H*log(10)) ./ log(Z);

    for k = 1:1000
        F  = 30*(delta/K).*log10(30*delta/K) - 1.2*H/K;
        dF = (30/K).*log10(30*delta/K) + (30/K)./log(10);

        delta_new = delta - F./dF;

        % Convergence check
        if max(abs(delta_new - delta)) < 1e-12
            delta = delta_new;
            converged = 1;
            fprintf('converged for %f\n',K)
            break;
        end
        delta = delta_new;
    end

    if converged == 0
        delta = NaN(size(delta));
    end


end


function fw = solve_for_fw(delta_ks)
    fw = 0.0604./(log10(22.*delta_ks).^2);
end


function [x, H, h,ub, A, xb] = shoaling_to_break_slope(H0, T, h0, m, g)


    npoints = 500;
    
    % x from offshore (0) to shoreline (h=0)
    x_max = h0 / m;
    x = linspace(0, x_max, npoints);
    
    h = h0 - m*x;   % depth along slope
    
    % dispersion and group velocity offshore
    k0 = waveNumber(T, h0,g);
    w = 2*pi/T;
    cg0 = groupVelocity(w, k0, h0);

    H = zeros(size(x));
    H(1) = H0;
    xb = NaN;

    ub = zeros(size(x));
    A = zeros(size(x));

    for i = 2:length(x)
        hi = h(i);
        if hi <= 0
            xb = x(i);
            break;
        end
    
        ki = waveNumber(T, hi,g);
        cg = groupVelocity(w, ki, hi);
    
        % shoaling first
        H(i) = H0 * sqrt(cg0 / cg);
    
        % now compute kinematics
        ub(i) = (pi*H(i)/T) / sinh(ki*hi);
        A(i)  = ub(i)/w;
    
        % breaking
        if H(i) >= 0.78 * hi
            xb = x(i);
            H(i) = 0.78 * hi;
            x = x(1:i);
            h = h(1:i);
            H = H(1:i);
            ub = ub(1:i);
            A  = A(1:i);
            break;
        end
    end

end


function k = waveNumber(T, h,g)
% Solve dispersion relation: w^2 = g k tanh(kh)
    w = 2*pi/T;
    k = w^2/g;   % deep water initial guess
    
    for n = 1:30
        f = g*k*tanh(k*h) - w^2;
        df = g*tanh(k*h) + g*k*h*(sech(k*h))^2;
        k = k - f/df;
    end
end


function cg = groupVelocity(w, k, h)
% Group velocity for linear waves
    cg = 0.5*(w/k)*(1 + (2*k*h)/sinh(2*k*h));
end

function Re_p = calc_particle_Reynolds(rho_s,rho,g,D50,nu)
    spec_weight = (rho_s/rho) - 1;
    Re_p = sqrt(spec_weight*g.*(D50.^3))./nu;


end

function critical_shields = calculate_Parker_shields_curve(Re_particle)
    critical_shields = 0.5.*(0.22.*(Re_particle.^(-0.6)) + 0.06.*(10.^(-7.7.*(Re_particle.^(-0.6)))));
end

function [s1,s2,s3,s4] = calc_suspension_threshold(rho_s,rho,g,D50,nu)

    P = 0.8;%1.2;
    kappa = 0.4;
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
    u_star = settling_velocity;%;./(P*kappa);
    suspension_threshold = (rho.*(u_star.^2))./((rho_s - rho)*g.*D50);

    % figure;
    % plot(D_star,W_star,'--')
    % hold on;
    % 
    % 
    % set(gca,'YScale','log');
    % set(gca,'XScale','log')


    bagnold_sus = 0.64.*((rho.*(settling_velocity.^2))./((rho_s - rho).*g.*D50)); % bagnold1966
    parker_sus = (rho.*(settling_velocity.^2))./((rho_s - rho)*g.*D50); % from online slides of parker (https://studylib.net/doc/9626909/1d-sediment-transport-morphodynamics-with)
    for i = 1:numel(D_star)
        REP = sqrt(D_star);
        if REP(i) > 1 && REP(i) < 32
            K = 4*(REP(i)^(-2/3));
        elseif REP(i) > 32
            K = 0.4;
        else
            K = NaN;
        end
        van_rijn_1984(i) = sqrt(K)*parker_sus(i);
    end
    % plot(sqrt(D_star),bagnold_sus,'--r')
    % plot(sqrt(D_star),parker_sus,'--g')
    % plot(sqrt(D_star),van_rijn_1984,'--b')
    % legend('Dietrich','Bagnold','Parker','Van Rijn')
    % 
    % xlabel('$D_*$','Interpreter','latex')
    % ylabel('$\tau$','Interpreter','latex')

    s1 = suspension_threshold;
    s2 = bagnold_sus;
    s3 = parker_sus;
    s4 = van_rijn_1984;

    % title('suspension thresholds')
end

function theta = calc_shields_number(rho_s,rho,depth,u,D50,g,fw)
    
    %theta = (rho*Manning*(abs(u)^2))./((rho_s - rho)*g.*D50.*(depth^(1/3)));
    kappa = 0.4;
    z0 = D50/30;
    u_star = (kappa*abs(u))./log((0.37*depth)./(z0));
    theta = (rho*fw*(u_star.^2))./((rho_s - rho)*g.*D50);
   
end

% 
% tau_b_sand =  0.5.*max(fw_sand).*rho.*ub.^2;
% 
% d_star = logspace(-3, 10, numel(tau_b_sand));
% theta_c = zeros(size(d_star));
% 
% % Shields curve approximation (Soulsby 1997)
% for i = 1:length(d_star)
%     ds = d_star(i);
% 
%     if ds < 4
%         theta_c(i) = 0.1 * ds^-1;
%     elseif ds < 10
%         theta_c(i) = 0.14 * ds^-0.64;
%     else
%         theta_c(i) = 0.055 * (1 + 0.3*ds)^-1;
%     end
% end
% 
% 
% d_star_sand = d50(1) * ((R*g)/(nu^2))^(1/3);
% d_star_sand = d_star_sand.*ones(size(tau_b_sand));
% 
% figure;
% plot(d_star,theta_c,'-k','LineWidth',3)
% hold on
% plot(d_star_sand,tau_b_sand,'-r')
% 
% function [tau_star, tau_c, tau_star_crit] = shields_critical_curve(tau, d,R,g)
% 
% 
%     % Dimensionless grain size
%     d_star = d * (( (rho * R * g) / (1e-6)^2 )^(1/3));  
% 
%     % Shields curve approximation (Soulsby 1997 / van Rijn)
%     if d_star < 4
%         tau_star_crit = 0.1 * d_star^-1;
%     elseif d_star < 10
%         tau_star_crit = 0.14 * d_star^-0.64;
%     else
%         tau_star_crit = 0.055 * (1 + 0.3*d_star)^-1;
%     end
% 
%     % Shields parameter
%     tau_star = tau / (rho * R * g * d);
% 
%     % critical shear stress
%     tau_c = tau_star_crit * rho * R * g * d;
% 
% end