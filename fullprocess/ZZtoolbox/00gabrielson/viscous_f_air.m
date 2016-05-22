function f_mu = viscous_f_air(freqs, radius, tempK, pressure)
% Viscous f function in air as a function of temperature [K],
% pressure [Pa], tube radius, [m] and frequency [Hz]
%
% USAGE: f_mu = viscous_f_air(freqs, radius, temperature, pressure);
%
% If pressure is not entered, 101325 Pa will be used; if temperature is
% not entered, 20 deg C will be used.
if nargin < 4; pressure = 101325; end
if nargin < 3; tempK = 273.1 + 20; end
%
d_mu = viscous_penetration_air(freqs, tempK, pressure);
argg = (j - 1)*radius./d_mu;
f_mu = 2*besselj(1, argg)./(argg.*besselj(0, argg));
