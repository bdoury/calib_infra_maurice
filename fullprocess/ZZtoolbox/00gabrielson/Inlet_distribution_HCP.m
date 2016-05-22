function [xx, yy] = Inlet_distribution_HCP(HCP_aperture, HCP_order)
%
%     X,Y coordinates for inlets of hexagonal close-packed pipe system
%
%       The arrangement of inlets for the hexagonal close-packed array
%   is specific.  The number of primary pipes and the number of tertiary
%   pipes must both be six.  To locate the centers of the individual hex
%   cells, the variable HCP_order must be specified.  The "order" is the
%   number of hex cells attached to each primary pipe (which is the same as
%   the number of secondary pipes).
%   
%       The HCP_aperture is the total distance across the inlet
%   distribution.  The variables, L_pri, L_sec, and L_ter are used in this
%   routine; however, they may be slightly longer than the pipe lengths
%   used in the acoustical model.  The calculation here is for the actual
%   inlet location, whereas the actual pipe lengths span from manifold to
%   manifold and will be slightly shorter.
%
%   USAGE: [xx, yy] = Inlet_distribution_HCP(HCP_aperture, HCP_order);
%
%   Examples (HCP's at I99)
%      For a 36-m aperture, 108-inlet HCP  
%        HCP_aperture = 36;  HCP_order = 3;  % which results in...
%           N_pri_pipes = 6;    N_sec_pipes = 3;   N_ter_pipes = 6;
%           L_pri = 10.28;  L_sec = 5.14;  L_ter = 2.57;  % lengths [m]
%
%      For a 18-m aperture,  72-inlet HCP  (HCP_order = 2):
%        HCP_aperture = 18;  HCP_order = 2;  % which results in...
%           N_pri_pipes = 6;    N_sec_pipes = 2;   N_ter_pipes = 6;
%           L_pri = 3.6;  L_sec = 3.6;  L_ter = 1.8;  % lengths [m]
%
N_pri_pipes = 6;  N_sec_pipes = HCP_order;  N_ter_pipes = 6;
L_basic = (HCP_aperture/2)/(2*HCP_order + 1);
L_pri = L_basic*2*(HCP_order - 1);
L_sec = 2*L_basic;  L_ter = L_basic;
%
d_pri_angle = 2*pi/N_pri_pipes;
pri_angles = (((1:N_pri_pipes) - 1)*d_pri_angle).';
%
%   Locations of hex-cell centers at vertices
xpv = (L_pri + L_sec)*cos(pri_angles);
ypv = (L_pri + L_sec)*sin(pri_angles);
%
inter_fraction = 1/N_sec_pipes;
%
nn = 1;
xpp = zeros(N_pri_pipes*(N_sec_pipes - 1), 1);
ypp = xpp;
%
for jj = 1:N_pri_pipes
   kk = jj + 1;
   if jj == N_pri_pipes; kk = 1; end
   slope = (ypv(kk) - ypv(jj))/(xpv(kk) - xpv(jj));
   dx =(xpv(kk) - xpv(jj))*inter_fraction;
   %   Locations of intermediate hex-cell centers
   for ii = 1:(N_sec_pipes - 1)
       xpp(nn) = xpv(jj) + dx*ii;
       ypp(nn) = ypv(jj) + dx*ii*slope;
       nn = nn + 1;
   end
end
%
%   Locations of inlets for vertex hex-cells
d_ter_angle = 2*pi/N_ter_pipes;
ter_angles = (((1:N_ter_pipes) - 1)*d_ter_angle).';
del_xc = L_ter*cos(ter_angles);
del_yc = L_ter*sin(ter_angles);
%
xvi = zeros(N_pri_pipes*(N_ter_pipes - 1), 1);
yvi = xvi; 
%
nn = 0;
%
for jj = 1:N_pri_pipes
    mm = 1:N_ter_pipes;  mx = mod(jj + 3, 6) + 1;
    mm(mx) = [];
    xvi((1:5) + nn) = xpv(jj) + del_xc(mm);
    yvi((1:5) + nn) = ypv(jj) + del_yc(mm);
    nn = nn + N_pri_pipes - 1;
end
%
xvii = zeros(N_pri_pipes*(N_sec_pipes - 1)*(N_ter_pipes - 1), 1);
yvii = xvii; 
%
nn = 0;  pp = 1;
%
for jj = 1:N_pri_pipes
    mm = 1:N_ter_pipes;  mx = mod(jj + 4, 6) + 1;
    mm(mx) = [];
    for kk = 1:(N_sec_pipes - 1)
        xvii((1:5) + nn) = xpp(pp) + del_xc(mm);
        yvii((1:5) + nn) = ypp(pp) + del_yc(mm);
        nn = nn + N_pri_pipes - 1;  pp = pp + 1;
    end
end
%
%   Collect all hex-cell centers in xpc, ypc
xpc = [xpv; xpp];  ypc = [ypv; ypp];
%
%   Collect all inlet positions in xp, yp
xx = [xpc; xvi; xvii];  yy = [ypc; yvi; yvii];
%
% % %  Diagnostic plots
% % figure(1)
% % hhh = plot(xpv, ypv, 'ro', xpp, ypp, 'bx',...
% %     xvi, yvi, 'go', xvii, yvii, 'k+');
% % %
% % figure(2)
% % hhh = plot(xx, yy, 'ro');
% % %