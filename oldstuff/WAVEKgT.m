function k = wavekgt( f, H, g, T, R, tol )
%
%    k = wavekgt( f, H, g, T, R, tol )
%
%    Returns a column vector of wavenumbers [k] (1/m) at frequencies [f] (Hz)
%    using the linear wave dispersion relation for water depth H (m).  tol
%    is an optional error bound on the predicted f's (the default is 1e-4). 
%    T is surface tension in N/m = dynes/cm * .001: default is 0.074
%    R = water density in Kgm/m^3: default is

%    Note: [f] > 0 must increase monotonically from small to large.  This
%    routine is optimized for speed and works best when length(f) is large.

                             % some initial housekeeping

if nargin<6,  tol = 1e-4;  end         % default tolerance
% g = 9.80171; 
c = 2*pi*sqrt(H/g);    
f = c * f(:);  Nf = length(f);         % non-dimensionalize f
B = T/(R*g*H*H);                       % non-dimensionalize T
k = zeros(Nf,1);


    k(1) = f(1)^2;                     % use deep water limit as initial guess
    for n=1:Nf
    dk=1;
    if n>1,
    k(n)=k(n-1);
    end

    while ( abs(dk) ) > tol         % Newton-Raphson iteration loop
        t = 1;
        if k(n)<20,
        t = tanh( k(n) ) ;
        end
        dk = -(f(n).^2 - k(n).*(1+B.*k(n).^2)*t)...
             ./ ( 3*B.*k(n).^2*t+t  + k(n)*(1+B.*k(n).^2)*( 1 - t.^2 ) ) ;
        k(n) = k(n) - dk ;
    end
    k(n) = abs( k(n) ) ;               % f(k) = f(-k), so k>0 == k<0 roots
end

N = min( find( k>50 ) ) ;              % inform user if err > tol
if isempty(N),  N = Nf;  end
err = abs( f(2:N) - sqrt( k(2:N).*tanh(k(2:N)) ) )./f(2:N) ;
%if max( err ) > tol,  fprintf('\n WAVEK: error exceeds %g \n', tol),  end

k = k / H ;                            % give k dimensions of 1/meter
