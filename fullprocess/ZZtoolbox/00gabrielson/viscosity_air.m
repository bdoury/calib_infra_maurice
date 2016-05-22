function mu = viscosity_air(tempK)
% Viscosity of air [Pa s] as a function of temperature [K];
% Sutherland's approximation.
%
% USAGE: mu = viscosity_air(temperature);
%
T_0 = 273.1; S = 110.6; % [K]
mu_0 = 17.16e-6; % [Pa sec]
alpha = 1.5;
%
mu = mu_0*((tempK/T_0).^alpha).*(T_0 + S)./(tempK + S);
