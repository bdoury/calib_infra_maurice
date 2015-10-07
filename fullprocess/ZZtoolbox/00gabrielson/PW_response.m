function p_inside = PW_response(freqs, sound_speed, arrival_elev,...
    arrival_azim, xx, yy)
%
%       Plane wave response for multi-inlet pipe model
%
%       For frequencies in 'freqs' vector; at Tkelvin temperature; and
%   for 'arrival_elev' in radians.  If 'arrival_azim' has multiple
%   values, then the response is averaged over those azimuthal angles.
%
%   USAGE:  p_inside = PW_rosette(freqs, Tkelvin, sound_speed, ...
%                          arrival_elev, xx, yy);
%
%       xx and yy are the X- and Y-coordinate vectors that specify the
%   positions of the inlets.  While this routine can be used to compute the
%   plane-wave response for an arbitrary distribution of inlets (or
%   sensors, for that matter), it will often be used in conjunction with an
%   acoustical model for the "internal" response of a pipe system.  In that
%   case, the assumption is that all of the inlets have the same internal
%   response.  If this is not true, the inlet locations should be grouped
%   by internal response and this routine should be run for each group.
% 
Nfreqs  = length(freqs);
Nangles = length(arrival_azim);
Nelements = length(xx);
%
%  Pre-allocate matrix for manifold pressure
p_inside_mat = zeros(Nfreqs, Nangles);
%
%  Loop over azimuth angles
for mm = 1:Nangles
    %  Loop over frequencies
    for nn = 1:Nfreqs
        omega = 2*pi*freqs(nn);
        k_horiz = omega*cos(arrival_elev)/sound_speed;
%         p_acs = exp(1j*k_horiz*xx);
        p_acs = exp(1j*k_horiz*(cos(arrival_azim(mm))*xx + ...
            sin(arrival_azim(mm))*yy));
        p_inside_mat(nn, mm) = sum(p_acs, 1)/Nelements;
    end
%
end
%  Average response over azimuth
p_inside = sum(p_inside_mat, 2)/Nangles;