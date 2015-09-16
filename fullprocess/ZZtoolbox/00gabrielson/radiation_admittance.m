function Y_rad = radiation_admittance(density, sound_speed, radius, freq);
% The function, 'radiation_admittance', computes the radiation admittance
% for an unbaffled circular exit with 'radius' assuming a medium with
% 'density' and 'sound_speed'. The parallel mass/resistance model is
% used, which is asymptotically correct at low frequency and is
% well-behaved at high frequency.
%
% USAGE: Y_rad = radiation_admittance(rho, c, a, f);
%
area = pi*radius^2;
Y_rad = 1/(1.44*density*sound_speed/area) + ...
1./(1j*0.6*2*pi*freq*density*radius/area);
