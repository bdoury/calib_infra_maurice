function [SkLfft, mseSk, MSC, ciMSC] = ...
    spectralmatrixestimation(...
    y, ...
    beta, ...
    Lfft,...
    w)
%=====================================================
% Spectral matrix estimation
% (smoothed periodogram with Bartlett window)
%
% synopsis :
% [SkLfft, mseSk, MSC, ciMSC] = ...
%       spectralmatrixestimation(y, mode, alpha, nfft)
%
%=====================================================
% y:  signal array of dimension N x d
%     N = sample number
%     d = dimension of the random process
% alpha: W = N^beta
%    where W is the window length
%    default value : alpha = 0.2
% nfft: FFT length
%======================
% Output:
%     SkLfft: Spectral matrix of dimension MxM for the
%           Lfft normalized frequencies (0:Lfft-1)/Lfft
%     mseSk: mean square error on SkLfft
%     MSC: magnitude squarred coherence
%     ciMSC: confidence interval on MSC
%=====================================================
[N,d]     = size(y);
if nargin==1
    beta = 0.2;
    Lfft  = N;
    w     = 'hamming';
    a     = 0.05;
end
if nargin == 2,
    Lfft = N;
    w    = 'hamming';
    a     = 0.05;
end
if nargin == 3,
    w    = 'hamming';
    a     = 0.05;
end
if nargin == 4,
    a     = 0.05;
end
y         = y-ones(N,1)*mean(y);
M         = ceil(N^(beta));
%=================================
% select window for smoothing
% must be with sum 1
cdwindow  = sprintf('WNa=window(@%s,2*M+1);',w);
eval(cdwindow);
% must be with sum 1
WN        = WNa/sum(WNa);
%=== Pk: periodogram
Pk        = zeros(d,d,N);
Sk        = zeros(d,d,N);
Yk        = fft(y,N)/sqrt(N);
for i1 = 1:d
    for i2 = 1:d
        Pk(i1,i2,:) = Yk(:,i1) .* conj(Yk(:,i2));
    end,
end
%=========== smoothing on the entire perodogram ===========
%=== Sk: smoothed periodogram
% warning: filtfilt avoids the FIR-delay of M
% but induces a stronger smoothing using twice filtering
for i1=1:d
    for i2=1:d
        auxI        = squeeze(Pk(i1,i2,:));
%         auxIe       = [auxI(1:N);auxI;auxI(2:end)];
%         auxSk       = filtfilt(WN,1,auxIe);
%         Sk(i1,i2,:) = auxSk(N+1:2*N);
        Sk(i1,i2,:) = conv(auxI,WN,'same');
    end
end
%===== selecting the uniform grid of length LFFT ==========
% if LFFT is not a divisor of N, the grid is approximative
Ngridf    = fix(N*(0:Lfft-1)/Lfft)+1;
SkLfft    = Sk(:,:,Ngridf);
combi     = d*(d-1)/2;
%==== MSE: mean square error on the smoothed periodogram
% theoretically the sum of variance plus squarred bias
mseSk     = zeros(d,d,Lfft);
for ik=M:Lfft-M-1
    G_ik = squeeze(Sk(:,:,ik-M+1:ik+M+1));
    [smoothEG, smoothvarG, smoothmse] = ...
        smooth2order(G_ik,WN);
    mseSk(:,:,ik) =  smoothmse;
end
%===== computation of the cross-spectral correlation
CSF       = zeros(combi,Lfft);
cpd       = 0;
for id1=1:d-1
    D1    = squeeze(sqrt(SkLfft(id1,id1,:)))+eps;
    A1    = max(D1);
    for id2=id1+1:d
        cpd        = cpd+1;
        D2         = squeeze(sqrt(SkLfft(id2,id2,:)))+eps;
        A2         = max(D2);
        D12        = squeeze((SkLfft(id1,id2,:)));
        CSF(cpd,:) = D12 ./ (D1 .* D2);
        CSF(cpd,union(find(D1<A1/10000),find(D2<A2/10000)))=0;
    end
end
%==== MSC for Magnitude Squarred Coherence
MSC    = abs(CSF) .^2;
%==== confidence interval on MSC
if nargout == 4;
    ciMSC  = zeros(combi,2,Lfft);
    for ic = 1:combi
        for ik=M:Lfft-M-1
            ciMSC(ic,1,ik) = invcumulFunctionMSC(a/2,MSC(ic,ik),2*M+1);
            ciMSC(ic,2,ik) = invcumulFunctionMSC(1-a/2,MSC(ic,ik),2*M+1);
        end
    end
end
%=========================================================================
%=========================================================================
function [smoothEG, smoothvarG, smoothmse] = ...
    smooth2order(Gs,WN)
%=========================================================================
twoMplus1  = length(WN);
M          = (twoMplus1-1)/2;
d          = size(Gs,1);
smoothEG   = zeros(d);
smoothvarG = zeros(d);

for im = 1:twoMplus1
    G           = squeeze(Gs(:,:,im));
    diagG       = diag(G);
    smoothEG    = smoothEG + G*WN(im);
    smoothvarG  = smoothvarG + diagG*diagG'*WN(im)*WN(im);
end
bias       = squeeze(Gs(:,:,M+1))-smoothEG;
smoothmse  = smoothvarG + bias*bias';
%=========================================================================
%=========================================================================
function c = invcumulFunctionMSC(p,C0,Nd,precision)
%=========================================================================
% Inverse cumulative function of the MSC estimates
% derived from the spectral estimation
%            |S_xy|^2
%    MSC = ------------
%           S_xx S_yy
%=========================================================================
% Synopsis
%          c = invcumulFunctionMSC(p,C0,Nd,precision)
% where
%   p
%   C0 : true coherence value in (0,1)
%   Nd : size of the smoothing window or averaging window
%
%=========================================================================
% outputs
%   P = Prob(\hat MSC =< E; N, MSC)
%   E 
%   CI : confidence interval at 100a%
%
% Rk: the resolution is obtained by dichotomy, whiic is easy
% because the cumulative function is by def monotonic.
%=========================================================================
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
%=========================================================================
%=========================================================================
function dP = cumulFunctionMSC(E,C0,Nd,p)
%=========================================================================
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
%=========================================================================

