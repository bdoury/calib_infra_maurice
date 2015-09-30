function [xx, yy] = Inlet_distribution_T_rosette(N_pri_pipes, L_pri_pipes,...
    N_sec_pipes, L_sec_pipes)
%
%       X,Y coordinates for inlets of a "T" rosette:  use this routine
%    with the acoustical model for an HCP.  This routine places two inlets
%    per secondary pipe in such a way that the inlets around each secondary
%    manifold are spaced with equal angles.
%
%  USAGE: [xx, yy] = Inlet_distribution_T_rosette(N_pri_pipes,...
%                      L_pri_pipes, N_sec_pipes, L_sec_pipes);
%
%
%  Angular orientation of secondary pipes
d_phi = 360/(N_sec_pipes);
d_phi_4 = d_phi/4;
pipe_phi_sec = (((1:N_sec_pipes) - 1)*d_phi + d_phi_4)*pi/180;
pipe_phi_sec = [pipe_phi_sec,...
    (((1:N_sec_pipes) - 1)*d_phi - d_phi_4)*pi/180];
%
%  Angular orientation of primary pipes
d_phi = 360/N_pri_pipes;
pipe_phi = (0:d_phi:359)*pi/180;
%
%  x, y coordinates of secondary manifolds relative to system center
xxp = L_pri_pipes*cos(pipe_phi);
yyp = L_pri_pipes*sin(pipe_phi);
%
%     Set up x,y coordinates for all inlets by adding primary offsets to
%  the secondary coordinates
Nelements = N_pri_pipes*N_sec_pipes;
xx = zeros(Nelements, 1);
yy = xx;  m1 = 1;
for ii = 1:N_pri_pipes
    for jj = 1:(2*N_sec_pipes)
    %  x,y coordinates relative to secondary manifold
    xxs = L_sec_pipes*cos(pipe_phi_sec(jj) + pipe_phi(ii));
    yys = L_sec_pipes*sin(pipe_phi_sec(jj) + pipe_phi(ii));
    xx(m1) = xxp(ii) + xxs;
    yy(m1) = yyp(ii) + yys;
    m1 = m1 + 1;
    end
end
%
m1 = 1;
xxt = zeros(2*N_sec_pipes, 1);  yyt = xxt;
%
for ii = 1:N_pri_pipes
    for jj = 1:N_sec_pipes
        xxt(jj) = xx(m1);  yyt(jj) = yy(m1);
        m1 = m1 + 1;
    end
    for jj = (N_sec_pipes + 1):(2*N_sec_pipes)
        xxt(jj) = xx(m1);  yyt(jj) = yy(m1);
        m1 = m1 + 1;
    end
%     for jj = 1:N_sec_pipes
%         xmid = (xxt(jj) + xxt(jj + N_sec_pipes))/2;
%         ymid = (yyt(jj) + yyt(jj + N_sec_pipes))/2;
%         hh = plot([xxp(ii), xmid], [yyp(ii), ymid]);
%         set(hh, 'linewidth', 2)
%         hh = plot([xxt(jj), xxt(jj + N_sec_pipes)],...
%             [yyt(jj), yyt(jj + N_sec_pipes)]);
%         set(hh, 'linewidth', 2)
%     end
end