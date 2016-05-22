function Zc = Z_capillary_inlet(freq, radius, Lc, tempK, pressure)
% Acoustical impedance of a capillary inlet with 'radius' and length,
% Lc as a function of frequency [Hz], temperature [K],
% and pressure [Pa]. Function returns a complex value.
% Frequency must be a scalar.
%
% USAGE: Zc = Z_capillary_inlet(freq, radius, Lc, temperature, pressure);
% (temperature and pressure are optional; if not entered, pressure =
% 101325 Pa and temperature = 273.1 + 20 K)
%
MW_air = 28.97; % [kg/kmole]
RR = 8314.3; % [J/kmole K]
if nargin < 5; pressure = 101325; end
if nargin < 4; tempK = 273.1 + 20; end
%
area = pi*radius^2;
Lcc = Lc + 2*0.6*radius; % acoustical mass correction, both ends
%
rho = pressure*MW_air/(RR*tempK);
mu = viscosity_air(tempK);
%
realZ = 8*mu*Lc/(pi*radius^4);
% % realZ = 2*realZ; % fudge
imagZ = 2*pi*freq*rho*Lcc/area; % use corrected length
%
Zc = realZ + 1j*imagZ;
