function T_rad = radiation_transferM(density, sound_speed, radius, freq);
% The function, 'radiation_transferM', computes the transfer matrix
% for radiation impedance for an unbaffled circular exit with 'radius'.
% A medium with 'density' and 'sound_speed' is assumed.
%
% USAGE: T_rad = radiation_transferM(rho, c, a, f);
%
Y_rad = radiation_admittance(density, sound_speed, radius, freq);
T_rad = [1 1/Y_rad; 0 1];
