function [p_inside, xx,yy] = PW_70_meter_I57(freqs, Tcelcius, arrival_theta,...
    N_pri_pipes, R_pri_pipes, N_sec_pipes, R_sec_pipes,...
    d_fuzz_phi)
%
% Plane wave response for rosette pipe model: 70-m rosette
% For frequencies in 'freqs' vector; at Tcelcius temperature; and for
% 'arrival_theta' in radians. 'arrival_theta' can be a vector
% containing several arrival angles.
%
% USAGE: p_inside = PW_70_meter_I57(freqs, Tcelcius, arrival_theta,...
% N_pri_pipes, R_pri_pipes, N_sec_pipes, R_sec_pipes,...
% d_fuzz_phi);
%
%
tempK = Tcelcius + 273.1; % ambient temperature [K]
cc    = 20.05*sqrt(tempK); % corresponding sound speed [m/s]
if nargin < 8
    d_fuzz_phi = 0;
end
%
% Angular orientation of secondary pipes
d_phi = 360/N_sec_pipes;
pipe_phi = (0:d_phi:359)*pi/180;
% x,y coordinates relative to secondary manifold
xxs = R_sec_pipes*cos(pipe_phi);
yys = R_sec_pipes*sin(pipe_phi);
%
% Angular orientation of primary pipes
% d_fuzz_phi allows perturbing radial locations of secondary manifolds
d_phi = 360/N_pri_pipes;
pipe_phi = (0:d_phi:359)*pi/180;
%
phi_fuzz = 2*d_fuzz_phi*(rand(size(pipe_phi)) - 0.5);
pipe_phi = pipe_phi + phi_fuzz*pi/180;
% x, y coordinates of secondary manifolds relative to system center
xxp = R_pri_pipes*cos(pipe_phi);
yyp = R_pri_pipes*sin(pipe_phi);
%
% Set up x,y coordinates for all inlets by adding primary offsets to the
% secondary coordinates
Nelements = N_pri_pipes*N_sec_pipes;
xx = zeros(Nelements, 1);
yy = xx; 
m1 = 1; 
m2 = N_sec_pipes;
for ii = 1:N_pri_pipes
    xx(m1:m2) = xxp(ii) + xxs;
    yy(m1:m2) = yyp(ii) + yys;
    m1 = m1 + N_sec_pipes; 
    m2 = m2 + N_sec_pipes;
end
%
Nfreqs = length(freqs);
Nangles = length(arrival_theta);
%
% Pre-allocate matrix for manifold pressure
p_inside = zeros(Nfreqs, Nangles);
%
% Loop over arrival angles
for mm = 1:Nangles
    % Loop over frequencies
    for nn = 1:Nfreqs
        omega = 2*pi*freqs(nn);
        k_horiz = omega*cos(arrival_theta(mm))/cc;
        p_acs = exp(1j*k_horiz*xx);
        p_inside(nn, mm) = sum(sum(p_acs, 2), 1)/Nelements;
    end
    %
end
