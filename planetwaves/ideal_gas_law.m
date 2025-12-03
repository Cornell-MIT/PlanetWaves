function rhoa = ideal_gas_law(Planet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the atmospheric density (kg/m3) using the ideal gas law
%           Ideal gas law: PV = nRT
% INPUT:
%   Planet : class containing information on 
%               surface_press    = surface pressure [Pa]
%               kgmolwt          = gram molecular weight of gas [Kgm/mol]
%               surface_temp     = surface temperature [K]
% OUTPUT:
%   rhoa   : atmospheric density [kg/m3]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RRR = 8.314;                                                               % Universal gas constant [J/K/mol]
    rhoa = Planet.surface_press*Planet.kgmolwt/(RRR*Planet.surface_temp);      % air density [kg/m3]
    %fprintf('rhoa in ideal gas law function: %f kg/m3\n', rhoa)
end