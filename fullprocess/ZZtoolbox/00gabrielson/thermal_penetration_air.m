function d_k = thermal_penetration_air(freqs, tempK, pressure)
% Thermal penetration depth [m] in air as a function of temperature [K],
% pressure [Pa], and frequency
%
% USAGE: d_k = thermal_penetration_air(freqs, temperature, pressure);
%
% If pressure is not entered, 101325 Pa will be used; if temperature is
% not entered, 20 deg C will be used.
if nargin < 3; pressure = 101325; end
if nargin < 2; tempK = 20 + 273.1; end
MW_air = 28.97; % [kg/kmole]
RR = 8314.3; % [J/kmole K]
cp = 1006; % [J/kg K]
kth = thermal_conductivity_air(tempK);
density = pressure*MW_air/(RR*tempK);
%
d_k = sqrt(2*kth./(density*2*pi*freqs*cp));
