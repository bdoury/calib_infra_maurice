function c = invcumulFunctionMSC(p,C0,Nd,precision)
%=============================================================
% Inverse cumulative function of the MSC estimates
% derived from the spectral estimation
%            |S_xy|^2
%    MSC = ------------
%           S_xx S_yy
%=============================================================
% Synopsis
%          c = invcumulFunctionMSC(p,C0,Nd,precision)
% where
%   p
%   C0 : true coherence value in (0,1)
%   Nd : size of the smoothing window or averaging window
%   precision: default value 1e-15
%=============================================================
% outputs
%   P = Prob(\hat MSC =< E; N, MSC)
%   E 
%   CI : confidence interval at 100a%
%
% Rk: the resolution is obtained by dichotomy, whiic is easy
% because the cumulative function is by def monotonic.
%=============================================================
if nargin < 4
    precision = 1e-15;
end
lastpos = 1;
lastneg = 0;
c     = 1;
while abs(lastpos-lastneg)>precision
    dy = cumulFunctionMSC(c,C0,Nd,p);
    if dy < 0
        lastneg = c;
        c = (lastneg+lastpos)/2;
    else
        lastpos = c;
        c = (lastneg+lastpos)/2;
    end
end
%=============================================================
function dP = cumulFunctionMSC(E,C0,Nd,p)
z = E*C0;
sum_F = 1;
coef_ell=1;
R = ((1-C0) ./ (1-z)) .^ Nd;
for ell=1:Nd-2
    Tk=1;
    sum_ell = 1;
    for k=1:ell
        Tk = (k-1-ell)*(k-Nd)* (Tk .*z) /k/k;
        sum_ell = sum_ell + Tk;
    end
    coef_ell = coef_ell .* ((1-E) ./ (1-z));
    sum_F = sum_F+coef_ell .* sum_ell;
end
dP = (R .* sum_F) .* E - p;
%=============================================================

