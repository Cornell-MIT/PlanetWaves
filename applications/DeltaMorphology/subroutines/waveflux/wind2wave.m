function [wave_height,wave_period,wave_length] = wind2wave(wind_speed,fetch,g)

% takes in properties of wind and computes steady state wave properties
% for volumetric flux calculations
    persistent warning_shown

    if isempty(warning_shown)
        warning('wind2wave: Using Fetch and Duration Limited Growth empirical functions, need to update to our model')
        warning_shown = true;
    end

    make_plot = 0;
    for i = 1:numel(wind_speed)
        u = wind_speed(i);
        for j = 1:numel(fetch)
            x = fetch(j);
            
            % dimensionless fetch
            X = (g*x)/(u^2); % eqn. 5.2 in Young+1991, Fetch and Duration Limited Growth

            % SMB empricial functions (eqn. 5.21-5.22 in Young+1991, Fetch and Duration Limited Growth)
            calc_H_dimless = @(X) (5e-3).*((tanh(0.0125.*(X.^0.42))).^2);
            calc_fp_dimless = @(X) (0.133.*((tanh(0.077.*(X.^0.25)))));
            
            
            Hs(i,j) = 4.*(((calc_H_dimless(X).^0.5).*(u.^2))./g); % eqn. 5.1
            fp = (calc_fp_dimless(X).*g)./u; % eqn. 5.2
            Tp(i,j) = 1./fp;

            Lp(i,j) = (g.*Tp(i,j).^2)./ (2*pi); % assume deepwater
    
        end
    end

    wave_height = Hs;
    wave_period = Tp;
    wave_length = Lp;
    
    if make_plot

          figure;
          tiledlayout('horizontal')
          nexttile
          plot(fetch./1000)
          xlabel('fetch km')
          nexttile
          plot(wave_height)
          xlabel('wave height')
          nexttile
          plot(wave_period)
          xlabel('wave period')


    end
end