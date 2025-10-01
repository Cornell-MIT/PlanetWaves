function [sand_river,gravel_river] = riverine_flux(rho_s,rho,nu,g,gravel_D50,fines_D50,width,slope)
% riverine sediment flux based on Birch+2023 
% INPUTS:
%   rho_s  = sediment density [kg/m3]
%   rho    = liquid density [kg/m3]

    if isnan(width)
        mode = 'slope-derivation';
        %disp('Computing non-dimensional flow parameters using slope-estimate')
    end
    if isnan(slope)
        mode = 'width-derivation';
        %disp('Computing non-dimensional flow parameters using width-estimate')
    end

    
    R = rho_s/rho - 1;
    Re_p = fines_D50*sqrt((R*g*fines_D50)/nu);
    

    if strcmp(mode,'width-derivation')
        bedload_constants = constants_bedload(R);
        gravel_river = bedload(g,width,gravel_D50,bedload_constants);

        susload_constants = constants_suspended_load(R,Re_p);
        sand_river = sus_load(g,width,fines_D50,R,nu,susload_constants);

    else
        error('slope derivation not implemented')
    end

% -- BEDLOAD -------------------------------------------------------------------------- %
    function gravel_river = bedload(g,width,D50,constants)

        alpha_b = constants.alpha_b;
        alpha_h = constants.alpha_h;
        alpha_s = constants.alpha_s;
        alpha_y = constants.alpha_y;

        n_b = constants.n_b;
        n_h = constants.n_h;
        n_s = constants.n_s;
        n_y = constants.n_y;

        Q = (((1/alpha_b)*(width*(g^(0.5*n_b + 0.2))*(D50^(2.5*n_b)))))^(1/(n_b + 0.4)); % eqn 4aa
        H = alpha_h*(Q^(n_h + 0.4))*(g^(-(0.5*n_h + 0.2)))*(D50^(-(2.5*n_h))); % eqn 10a
        B = alpha_b*((Q^(n_b + 0.4))*(g^(-(0.5*n_b + 0.2)))*(D50^(-2.5*n_b))); % eqn 10b
        S = alpha_s*((Q^n_s)*(g^(-0.5*n_s))*(D50^(-2.5*n_s))); % eqn 10c
        Qs = alpha_y*((Q^n_y)*(g^(0.5*(1-n_y)))*(D50^(2.5*(1-n_y)))); % eqn 10d

        % returns the values with SI units (not the non dimensional forms)
        gravel_river.B = B;
        gravel_river.H = H;
        gravel_river.S = S;
        gravel_river.Q = Q;
        gravel_river.Qs = Qs;

    end

    function sand_river = sus_load(g,width,D50,R,nu,constants)
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

        Q = (((width/alpha_b))*(g^(0.5*(n_b - m_b) + 0.2))*(D50^(2.5*n_b - 1.5*m_b))*(R^(-0.5*m_b))*(nu^m_b))^(1/(n_b + 0.4)); % eqn 5a
        
        H = alpha_h*(Q^(n_h + 0.4))*(g^(0.5*(m_h - n_h) - 0.2))*(D50^(1.5*m_h - 2.5*n_h))*(R^(0.5*m_h))*(nu^(-m_h)); % eqn 11a
        B = alpha_b*(Q^(n_b + 0.4))*(g^(-(0.5*n_b + 0.2) + 0.5*m_b))*(D50^(-2.5*n_b + 1.5*m_b))*(nu^-m_b)*(R^(0.5*m_b)); % eqn 11b
        S = alpha_s*(Q^n_s)*(g^(0.5*(m_s - n_s)))*(D50^(-2.5*n_s + 1.5*m_s))*(nu^(-m_s))*(R^(0.5*m_s)); % eqn 11c
        Qs = (3.6e-5)*(Q^n_y)*(g^(0.5*(1-n_y) -0.06))*(D50^(2.5*(1-n_y) - 0.18))*(nu^0.12)*(R^-0.06);

        sand_river.B = B;
        sand_river.H = H;
        sand_river.S = S;
        sand_river.Q = Q;
        sand_river.Qs = Qs;
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

