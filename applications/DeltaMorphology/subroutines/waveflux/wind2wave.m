function [wave_height,wave_period,wave_length] = wind2wave(wind_speed,fetch,g)

% takes in properties of wind and computes steady state wave properties
% for volumetric flux calculations

    persistent warning_shown

    if isempty(warning_shown)
        warning('wind2wave: Using SMB/Young(1991) Fetch-Limited Growth with JONSWAP comparison.')
        warning_shown = true;
    end

    make_plot = 1;
    compare2jonswap = 1;

    wind_speed = wind_speed(:)';
    fetch = fetch(:)';

    wave_height = NaN(numel(wind_speed), numel(fetch));
    wave_period = NaN(numel(wind_speed), numel(fetch));
    wave_length = NaN(numel(wind_speed), numel(fetch));

    % SMB Young
    for i = 1:numel(wind_speed)
        u = wind_speed(i);
        X = (g*fetch)/(u^2); % dimensionless fetch

        % Empirical dimensionless functions
        Hs_dimless = (5e-3) .* (tanh(0.0125.*(X.^0.42))).^2;
        fp_dimless = 0.133 .* tanh(0.077.*(X.^0.25));

        % Convert to dimensional parameters
        Hs = 4 .* sqrt(Hs_dimless) .* (u.^2./g);
        fp = fp_dimless .* (g./u);
        Tp = 1./fp;
        Lp = (g.*Tp.^2)./(2*pi); % deep-water wavelength

        wave_height(i,:) = Hs;
        wave_period(i,:) = Tp;
        wave_length(i,:) = Lp;
    end

    % compare with other empircal relationship
    if compare2jonswap
        for i = 1:numel(wind_speed)
            u = wind_speed(i);
            X = (g*fetch)/(u^2);

            % JONSWAP/SPM empirical formulas 
            Hs_J = 0.283 .* (tanh(0.0125*X.^0.42).*tanh(0.0125*X.^0.75)) .* (u.^2/g);
            Tp_J = 7.54 .* tanh(0.077.*X.^0.25) .* (u./g);

            % Compute relative differences at middle fetch range
            mid_idx = round(numel(fetch)/2);
            H_diff = abs(wave_height(i,mid_idx) - Hs_J(mid_idx)) / Hs_J(mid_idx);
            T_diff = abs(wave_period(i,mid_idx) - Tp_J(mid_idx)) / Tp_J(mid_idx);

            if H_diff > 0.3
                warning('Hs differs from JONSWAP by >30%% for U=%.1f m/s (mid fetch ~%.0f m)', ...
                    u, fetch(mid_idx));
            end
            if T_diff > 0.2
                warning('Tp differs from JONSWAP by >20%% for U=%.1f m/s (mid fetch ~%.0f m)', ...
                    u, fetch(mid_idx));
            end

            if make_plot
                figure(1); 
                subplot(2,1,1); hold on;
                loglog(fetch/1000, wave_height(i,:), '-', 'LineWidth',1.6, 'DisplayName', sprintf('SMB %g m/s',u));
                loglog(fetch/1000, Hs_J, '--', 'LineWidth',1.6, 'DisplayName', sprintf('JONSWAP %g m/s',u));
                xlabel('Fetch (km)'); ylabel('Hs (m)'); grid on; title('Hs vs Fetch');

                subplot(2,1,2); hold on;
                loglog(fetch/1000, wave_period(i,:), '-', 'LineWidth',1.6, 'DisplayName', sprintf('SMB %g m/s',u));
                loglog(fetch/1000, Tp_J, '--', 'LineWidth',1.6, 'DisplayName', sprintf('JONSWAP %g m/s',u));
                xlabel('Fetch (km)'); ylabel('Tp (s)'); grid on; title('Tp vs Fetch');
            end
        end

        if make_plot
            subplot(2,1,1); legend show;
            subplot(2,1,2); legend show;
        end
    end


    
end