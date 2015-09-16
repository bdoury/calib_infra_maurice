function d_mu = viscous_penetration_air(freqs, tempK, pressure)
% Viscous penetration depth [m] in air as a function of temperature [K],
% pressure [Pa], and frequency
%
% USAGE: d_mu = viscous_penetration_air(freqs, temperature, pressure);
%
% If pressure is not entered, 101325 Pa will be used; if temperature is
% not entered, 20 deg C will be used.
if nargin < 3; pressure = 101325; end
if nargin < 2; tempK = 20 + 273.1; end
MW_air = 28.97; % [kg/kmole]
RR = 8314.3; % [J/kmole K]
%
mu = viscosity_air(tempK);
density = pressure*MW_air/(RR*tempK);
%
d_mu = sqrt(2*mu./(density*2*pi*freqs));
