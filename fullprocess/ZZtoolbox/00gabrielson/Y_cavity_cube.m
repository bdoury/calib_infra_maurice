function Ym = Y_cavity_cube(freq, length, tempK, pressure)
% Acoustical admittance of a cavity with 'radius' and 'length',
% as a function of frequency [Hz], temperature [K],
% and pressure [Pa]. Function returns a complex value.
% Frequency must be a scalar. Based on scaled results for sphere adapted
% for cylindrical geometry and valid for any frequency.
%
% USAGE: Ym = Y_cavity(freq, radius, length, temperature, pressure);
% (temperature and pressure are optional; if not entered, pressure =
% 101325 Pa and temperature = 273.1 + 20 K)
%
gamma = 1.40; % ratio of specific heats
if nargin < 5; pressure = 101325; end
if nargin < 4; tempK = 273.1 + 20; end
% Geometry: assumed cylindrical volume
end_area = length^2;
volume = length^3;
compliance = volume/(gamma*pressure); % adiabatic compliance
total_area = 6*end_area;
% Calculate admittance of manifold
YC = 2j*pi*freq*compliance; % adiabatic admittance
%
d_kth = thermal_penetration_air(freq, tempK, pressure);
V_S_cyl = volume/total_area;
radius_sphere = 3*V_S_cyl; % sphere with equivalent radius
arg = (radius_sphere/d_kth)*(1 + 1j);
FF = 1 - (3/arg)*coth(arg) + 3/(arg^2);
%
Ym = YC*(gamma - (gamma - 1)*FF); % corrected admittance
