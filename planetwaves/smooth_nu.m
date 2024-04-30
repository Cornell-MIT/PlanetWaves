function ustar = smooth_nu(U,z,nu)
% To calculate the friction velocity for smooth flow of any gas.
%
% function ustar = smooth(U,z)
%
% ustar = friction velocity [m/s]
% U     = wind speed [m/s] at height z
% z     = height of wind speed measurement [m]
% nu    = kinematic viscosity of gas [m^2/s]

    ustar = NaN(1,length(U));
    z0 = 0.001;

    for j = 1:length(U)
       m = 0;    
       for k = 1:6
           m = m+1;
           ust = 0.4*U(j)/(log(z/z0));
           z0 = (1/9)*nu/ust;
       end
       ustar(j) = ust;
    end
end