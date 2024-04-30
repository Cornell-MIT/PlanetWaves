function k = wavekgt( f, D, g, T, rho, tol )
%
%    k = wavekgt( f, H, g, T, rho, tol )
%
%    Returns a column vector of wavenumbers [k] (1/m) at frequencies [f] (Hz)
%    using the linear wave dispersion relation for water depth H (m).  tol
%    is an optional error bound on the predicted f's (the default is 1e-4).
%    T is surface tension in N/m = dynes/cm * .001: default is 0.074
%    rho = water density in Kgm/m^3: default is
%    Note: [f] > 0 must increase monotonically from small to large.  This
%    routine is optimized for speed and works best when length(f) is large.

    k = zeros(length(f),1);               % initialize vector of wavenumbers to be filled            
    if D > 0
        if nargin<6
            tol = 1e-4;                   % default tolerance
        end                                                                                   
        f = (2*pi*sqrt(D/g)) * f(:);      % non-dimensionalize angular frequency
        B = T/(rho*g*D*D);                % non-dimensionalize surface tension

        k(1) = f(1)^2;                    % use deep water limit as initial guess
        for n = 1:length(f)
            dk = 1;
            if n > 1
                k(n) = k(n-1);
            end
            maxWhile = 1e8;
            numWhile = 0;
            while ( abs(dk) ) > tol       % Newton-Raphson iteration loop
                t = 1;                    % tanh(x) -> 1 as x -> infinity 
                if k(n) < 20              % argument of tanh(x) is considered large enough to aprox to 1 when x > 20
                    t = tanh(k(n)) ;
                end
                dk = -(f(n).^2 - k(n).*(1+B.*k(n).^2)*t)...
                    ./ ( 3*B.*k(n).^2*t+t  + k(n)*(1+B.*k(n).^2)*( 1 - t.^2 ) ); 
                k(n) = k(n) - dk ;
                numWhile = numWhile + 1;
                if numWhile > maxWhile
                    warning(sprintf('stuck in while loop calculating wave number. Using less tuned value with error: %0.5f',dk))
                    break;
                end
            end
            k(n) = abs(k(n)) ;               
        end
    
        k = k/D ;                          % give k dimensions of 1/meter
    else
        k = NaN(length(f),1);
    end
end