function TT = Tmatrix_air(freq, radius, dx, tempK, pressure, correction)
% Transfer matrix for an incremental length, dx [m], of circular pipe
% in air as a function of frequency [Hz], radius [m], temperature [K],
% and pressure [Pa]. Function returns a 2x2 matrix. Frequency must be a
% scalar. 'correction' is an empirical correction factor for the
% effective damping. Normally, it would be set to one (and is set to one
% if the argument is omitted). Set to zero to disable the correction.
%
% USAGE: TT = Tmatrix_air(freq, radius, dx, kelvin_temp, pressure);
% (temperature and pressure are optional; if not entered, pressure =
% 101325 Pa and temperature = 20 deg C)
%
gamma = 1.40;
MW_air = 28.97; % [kg/kmole]
RR = 8314.3; % [J/kmole K]
if nargin < 6; correction = 1; end
if nargin < 5; pressure = 101325; end
if nargin < 4; tempK = 273.1 + 20; end
%
f_mu = viscous_f_air(freq, radius, tempK, pressure);
f_k = thermal_f_air(freq, radius, tempK, pressure);
area = pi*radius^2;
%
density = pressure*MW_air/(RR*tempK);
ZA = 2j*pi*freq*density*dx./(2*area*(1 - f_mu));
YB = 2j*pi*freq*dx*area*(1 + (gamma - 1)*f_k)/(gamma*pressure);
ZA = (correction + 1)*real(ZA) + 1j*imag(ZA); % empirical correction
YB = (correction + 1)*real(YB) + 1j*imag(YB); % empirical correction
%
Zrat = ZA.*YB;
TT = [1 + Zrat, ZA.*(2 + Zrat); YB, 1 + Zrat];

