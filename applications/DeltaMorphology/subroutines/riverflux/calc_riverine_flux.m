function [suspendedload_dominated,bedload_dominated] = calc_riverine_flux(planet,river)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% riverine sediment flux based on Birch+2023 based on river's width and slope
% INPUTS:
%   planet
%       rho_s  = sediment density [kg/m3]
%       rho    = liquid density [kg/m3]
%       g      = gravity [m/s^2]
%       nu     = viscosity [m^2/s]
%   river
%       width      = width of river cross-stream [m]
%       slope      = slope of river downstream [m/m]
% OUTPUTS:
%   bedload_dominated
%       D50    = grain size [m]
%       H      = river depth [m]
%       Q      = flow discharge [m3/s]
%       Qs     = sediment flux [m3/s]
%   suspendedload_dominated
%       D50    = grain size [m]
%       H      = river depth [m]
%       Q      = flow discharge [m3/s]
%       Qs     = sediment flux [m3/s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    R = planet.rho_s/planet.rho - 1;

    % calculate river discharge, sediment size, and sediment flux for bedload dominated
    bedload_dominated = calc_bedload_dominated(river.width,river.slope,planet.g,constants_bedload(R));
    
    Re_p = (bedload_dominated.D50*sqrt(R*planet.g*bedload_dominated.D50))/planet.nu;
    % calculate river discharge, sediment size, and sediment flux for suspended load dominated river
    suspendedload_dominated = calc_susload_dominated(river.width,river.slope,planet.g,planet.nu,R,constants_suspended_load(R,Re_p));

% -- BEDLOAD -------------------------------------------------------------------------- %
    function bedload_dominated = calc_bedload_dominated(B,S,g,constants)

        alpha_b = constants.alpha_b;
        alpha_h = constants.alpha_h;
        alpha_s = constants.alpha_s;
        alpha_y = constants.alpha_y;

        n_b = constants.n_b;
        n_h = constants.n_h;
        n_s = constants.n_s;
        n_y = constants.n_y;

        % returns the values with SI units (not the non dimensional forms)
        D50 = (B/alpha_b)*((alpha_s/S)^((n_b + 0.4)/n_s)); % Eqn. 4bb
        Q = ((S/alpha_s)^(1/n_s))*(g^0.5)*(D50^2.5); % Eqn. 4aa
        H = alpha_h*(Q^(0.4+n_h)*g^(-0.2 - 0.5*n_h)*(D50^(-2.5*n_h)));

        Q_nondim = Q/((g^0.5)*(D50^2.5));

        Qs = alpha_y*(Q_nondim^n_y);

        bedload_dominated.D50 = D50;
        bedload_dominated.H = H;
        bedload_dominated.Q = Q;
        bedload_dominated.Qs = Qs;
        

    end

    function susload_dominated = calc_susload_dominated(B,S,g,nu,R,constants)
        alpha_b = constants.alpha_b;
        alpha_h = constants.alpha_h;
        alpha_s = constants.alpha_s;
        alpha_y = constants.alpha_y;
        

        m_b = constants.m_b;
        m_h = constants.m_h;
        m_s = constants.m_s;

        n_b = constants.n_b;
        n_h = constants.n_h;
        n_s = constants.n_s;
        n_y = constants.n_y;

        exp1 = (n_b + 0.4)/n_s;
        exp2 = -0.5*n_b*(m_s/n_s) - 0.2*(m_s/n_s) + 0.5*m_b;
        exp3 = 1/(1.5*n_b*(m_s/n_s) + 0.6*(m_s/n_s) - 1.5*m_b - 1);
        D50 = ((alpha_b/B)*((S/alpha_s)^(exp1))*(((R*g)/(nu^2))^exp2))^exp3; % Eqn. 5b
        Q = ((B/alpha_b)*(g^(0.5*(n_b - m_b)+0.2))*(D50^(2.5*n_b - 1.5*m_b))*(R^(-0.5*m_b))*(nu^m_b))^(1/(n_b + 0.4)); % Eqn. 5a
        H = alpha_h*(Q^(n_h + 0.4)*(g^(0.5*(m_h - n_h)-0.2)*(D50^(1.5*m_h - 2.5*n_h))*(R^(0.5*m_h))*(nu^(-m_h)))); 
        Qs = alpha_y*(Q^n_y)*(g^((1-n_y)/2))*(D50^(2.5*(1-n_y)));


        susload_dominated.D50 = D50;
        susload_dominated.H = H;
        susload_dominated.Q = Q;
        susload_dominated.Qs = Qs;
    end

% -- POWER LAW CONSTANTS ---------------------------------------------------------------- %
    % SUSPENDED LOAD
    function susload = constants_suspended_load(R,Re_p)
        % Table S2
        susload.alpha_b = 0.95/(sqrt(R)); % eqn 3
        susload.n_b = 0.11; % Table S2
        susload.m_b = 0.1; % Table S2

        susload.alpha_h = 3.8; % eqn 3
        susload.n_h = -0.06; % Table S2
        susload.m_h = -0.11; % Table S2

        susload.alpha_s = 0.01*R; % eqn 3
        susload.n_s = -0.17; % Table S2
        susload.m_s = -0.05; % Table S2

        susload.alpha_y = (3.6e-5)*(Re_p^-0.12); % equation 3 (eqn 33 in supplement)
        susload.n_y = 1.1; % equation 31 in supplement
    
    end
    % BEDLOAD
    function bedload = constants_bedload(R) 
        
        bedload.alpha_b = 18/((1+R)*sqrt(R)); % eqn 2
        bedload.n_b = 0.06; % Table S2
        bedload.m_b = NaN; % Table S2

        bedload.alpha_h = 0.22*((1+R)^0.73); % eqn 2
        bedload.n_h = -0.02; % Table S2
        bedload.m_h = NaN; % Table S2

        bedload.alpha_s = (0.11*R)/((1+R)^0.73); % eqn 2
        bedload.n_s = -0.33; % Table S2
        bedload.m_s = NaN; % Table S2

        bedload.alpha_y = 0.01/(1+R); % eqn 2
        bedload.n_y = 0.53; % equation 23 in supplement

    end

end

