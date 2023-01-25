function ustar = smooth_nu(U,z,nu) 

% To calculate the friction velocity for smooth flow of any gas.
%
% function ustar = smooth(U,z) 
%
% ustar = friction velocity [m/s]
% U     = wind speed [m/s] at height z
% z     = height of wind speed measurement [m]
% nu    = kinematic viscosity of gas [m^2/s]
%


len=length(U);
z0=0.001;
z1=0.002;
% nu=0.156/10000;
del=0.00001;

for j=1:len
m=0;
u=U(j);
z1=z0+2*del;
%while abs(z0-z1) > del

for k=1:6
m=m+1;
z1=z0;
ust=0.4*u/(log(z/z0));
z0=0.132*nu/ust;
end
ustar(j)=ust;
end
