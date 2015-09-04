function P = cumulFunctionMSC(E,MSC,N)
%=============================================================
% Cumulative function of the MSC estimator
% derived from the spectral estimation
%            |S_xy|^2
%    MSC = ------------
%           S_xx S_yy
%=============================================================
% Synopsis
%          P = cumulFunctionMSC(MSC,N)
% where
%   MSC : true coherence value in (0,1)
%   N : size of the smoothing window
%
%=============================================================
% outputs
%   P = Prob(\hat MSC =< E; N, MSC) 
%=============================================================
LE = length(E);
z = E*MSC;
sum_F = 1;
coef_ell=1;
R               = ((1-MSC) ./ (1-z)) .^ N;
for ell=1:N-2
    Tk          = ones(LE,1);
    sum_ell     = 1;
    for k=1:ell
        Tk      = (k-1-ell)*(k-N)* (Tk .*z) /k/k;
        sum_ell = sum_ell + Tk;
    end
    coef_ell    = coef_ell .* ((1-E) ./ (1-z));
    sum_F       = sum_F+coef_ell .* sum_ell;
end
P               = (R .* sum_F) .* E;
%=============================================================
