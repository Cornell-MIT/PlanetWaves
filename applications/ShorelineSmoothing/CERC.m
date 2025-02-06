function Qs = CERC(H,T,theta)
   % FIND SEDIMENT FLUX FOR WAVES OF HEIGHT (H), PERIOD (T), WITH ATTACK ANGLE (theta)
   % INPUT:
   %    H = wave height [m]
   %    T = wave period [s]
   %    theta = wavefront angle - shoreface angle 
   % OUTPUT:
   %    Qs = Volumetric Sediment Flux

    Qs = ((H).^(12/5)).*((T).^(1/5)).*(cos(theta).^(6/5)).*sin(theta);
end
