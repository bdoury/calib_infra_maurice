function pdf = pdfMSC(E,MSC,N)
%=============================================================
% probability density of the MSC estimates
% derived from the spectral estimation
%            |S_xy|^2
%    MSC = ------------
%           S_xx S_yy
%=============================================================
% Synopsis
%          pdf = pdfMSC(E,MSC,N)
% where
%   E : range in (0,1)
%   MSC : true coherence value in (0,1)
%   N : size of the smoothing window
%
%=============================================================
% outputs:
%   pdf: probability density function of MSC
%=============================================================
LE       = length(E);
Emu      = E*MSC;
R1       = (((1-MSC) .*(1-E)) ./ ((1-Emu) .^2)) .^ N;
R2       = (1-Emu) ./ ((1-E) .^2);
Tk       = ones(LE,1);
sumT     = 1;
for k=1:N-1
    Tk   = (k-N)*(k-N)* (Tk .*Emu) /k/k;
    sumT = sumT + Tk;
end
pdf      = (N-1) * (R1 .* R2) .* sumT;
%=============================================================
