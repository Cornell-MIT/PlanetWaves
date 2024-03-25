function [rho_liq,nu_liq,surface_ten] = lakecomp(methane, temp)
% Find the material properties of Titan lake given various compositons
% (Steckloff et al., 2020)
addpath('TitanData/Compositions')
%% LOAD DATA:
% top row = methane alkane fraction (molar fraction of the alkanes (CH4+C2H6) that is methane)
% first column = temperature [K]
format long
density = load('density.csv'); % [kg/m3]
kinematic_visc = load('kinematic_visc.csv'); % [cm2/s]
%dynamic_visc = load('Data/compositions/dynamic_visc.csv'); % [Pa*s]
surf_ten = load('surface_tension.csv'); % [N/m]
methane_frac = load('CH4_molefraction.csv'); % methane
%ethane_frac = load('Data/compositions/C2H6_molefraction.csv'); % ethane
%nitrogen_frac = load('Data/compositions/N2_molefraction.csv'); % nitrogen

temp_location = find(methane_frac(:,1)==temp);               % Find row with temperature

%nitrogen = nitrogen_frac(temp_location, max(find(abs(nitrogen_frac(1,:) - methane)<0.001))); % Fraction of nitrogen
rho_liq = density(temp_location, max(find(abs(density(1,:) - methane)<0.001))); % liquid density [kg/m3]
nu_liq = kinematic_visc(temp_location, max(find(abs(kinematic_visc(1,:) - methane)<0.001)))/10000;  % liquid viscocity [m2/s]
surface_ten = surf_ten(temp_location, max(find(abs(surf_ten(1,:) - methane)<0.001))); % liquid surface tension [N/m]

end