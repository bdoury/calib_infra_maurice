function kth = thermal_conductivity_air(tempK)
% Thermal conductivity [W/m K] of air as a function of temperature [K];
% Sutherland's approximation.
%
% USAGE: kth = thermal_conductivity_air(temperature);
%
T_0 = 273.1; S = 194.4; % [K]
kth_0 = 0.02413; % [W/m K]
alpha = 1.5;
%
kth = kth_0*((tempK/T_0).^alpha).*(T_0 + S)./(tempK + S);
