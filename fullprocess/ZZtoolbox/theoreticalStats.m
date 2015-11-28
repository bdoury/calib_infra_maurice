function [STATSmodUUonRU, STATSmodURonRR, statMSC, ...
    STDphase_rd, statPhasepsd] = ...
    theoreticalStats(allT, Smatrix, N, alphaCI)
%=========================================================================
% Compute the probability density functions, the STDs 
% and the CI (at alphaCI) of:
%    - the ratio modules:
%
%             |Smatrix(1,2)|                  Smatrix(1,1)
%    Rinf =  ---------------   and   Rsup =  --------------
%              Smatrix(2,2)                  |Smatrix(2,1)|
%    - the phase
%                phase = arg Smatrix(1,2)
%    - the MSC 
%                        |Smatrix(1,2)|^2
%                MSC = -------------------------
%                       Smatrix(1,1)Smatrix(2,2)
%=========================================================================
% Used Matlab functions
%    INTEGRAL on R2013
%    QUADGK   on R2010
%    besseli
%=========================================================================
% Inputs :
%       - allT:
%             allT.TUUonUR: list of values range of SUU/abs(SUR)
%             allT.TURonRR: list of values range of abs(SUR)/SRR
%             allT.phase:   list of values of arg(S12)
%             allT.MSC: list of values of MSC
%       - Smatrix: is the 2 x 2 spectral matrix at one frequency bin
%                   | Smatrix(1,1) Smatrix(1,2) |
%                   | Smatrix(2,1) Smatrix(2,2) |
%        Smatrix(1,1) is spectrum on SUT, Smatrix(2,2) is spectrum on SREF
%       - N: window length for averaging FFT frames
%       - alphaCI: level of confidence, typical values 0.05, 0.15, 0.3
%  for example for alphaCI = 2*(1-normcdf(1)) about 0.31,then CI = 2*sigma
%=========================================================================
% STATSmodUUonRU structure
%              pdf: probability density function
%              median
%              mean
%              cumul: cumulative function
%              CI: confidence interval
% STATSmodURonRR structure
%              pdf: probability density function
%              median
%              mean
%              cumul: cumulative function
%              CI: confidence interval
% statMSC structure, 
%              pdf: probability density function
%              CI: confidence interval
% STDphase_rd : std on the pahse
% statPhasepsd : probability density function of the phase
%=========================================================================
% values of MSC range
valMSC              = allT.MSC;
%=======
rho                 = abs(Smatrix(1,2)/sqrt(Smatrix(1,1)*Smatrix(2,2)));
MSC_theo            = rho ^2;
lambda              = sqrt(Smatrix(1,1)/Smatrix(2,2));
% performs (1*2*...*(N-1))^(1/N)
gammaNm1P           = exp(sum(log(1:N-1))/N);
phase_HUminusHR_rad = -atan2(imag(Smatrix(1,2)),real(Smatrix(1,2)));
%======================================================================
%====================== phase =========================================
%======================================================================
STDphase_rd  = asin(sqrt((1-MSC_theo)/N/MSC_theo/2));
statPhasepsd = ...
    (1/sqrt(2*pi)/STDphase_rd)* ...
    exp(-(allT.phase_rd-phase_HUminusHR_rad) .^2/...
    (2 * STDphase_rd^2));
%======================================================================
%====================== first ratio ===================================
%======================================================================
TURonRR         = allT.TURonRR;
TURonRR         = TURonRR(:);
zeta      = (1-MSC_theo) ./ (2*rho*lambda * TURonRR)/gammaNm1P;
cst       = lambda/rho;
xi        = ((1+lambda*lambda*(TURonRR .* TURonRR))) ./(2*rho*lambda*TURonRR);
xim1      = xi-1;
LT        = length(TURonRR);
pURonRR   = zeros(LT,1);
for it=1:LT
    if exist('integral','file')
        pURonRR(it)   = cst * integral(@(x) myf(x,xim1(it), ...
            zeta(it),N),0,inf);
    else
        pURonRR(it)   = cst * quadgk(@(x) myf(x,xim1(it), ...
            zeta(it),N),0,inf);
    end
end
STATSmodURonRR.pdf = pURonRR;
valMSC             = valMSC(:);
statMSC.pdf        = pdfMSC(valMSC,MSC_theo,N);
statMSC.CI(1)      = invcumulFunctionMSC(alphaCI/2,MSC_theo,N);
statMSC.CI(2)      = invcumulFunctionMSC(1-alphaCI/2,MSC_theo,N);
%======================================================================
%===================== second ratio ===================================
%======================================================================
TUUonUR   = allT.TUUonUR;
TUUonUR   = TUUonUR(:);
Tper      = 1 ./ TUUonUR;
lambda    = sqrt(Smatrix(2,2)/Smatrix(1,1));
zeta      = (1-MSC_theo) ./ (2*rho*lambda * Tper)/gammaNm1P;
cst       = lambda/rho;
xi        = ((1+lambda*lambda*(Tper .* Tper))) ./ (2*rho*lambda*Tper);
xim1      = xi-1;
LTper     = length(Tper);
pUUonRU   = zeros(LTper,1);
for it=1:LTper
    if exist('integral','file')
        pUUonRU(it)   = cst * integral(@(x) myf(x,xim1(it), ...
            zeta(it),N),0,inf);
    else
        pUUonRU(it)   = cst * quadgk(@(x) myf(x,xim1(it), ...
            zeta(it),N),0,inf);
    end
end
STATSmodUUonRU.pdf = pUUonRU .* (Tper .* Tper) ;
cumul12on22 = cumsum(STATSmodURonRR.pdf) * (TURonRR(2)-TURonRR(1));
STATSmodURonRR.mean = ...
    sum(STATSmodURonRR.pdf .* TURonRR) * (TURonRR(2)-TURonRR(1));
STATSmodURonRR.median = TURonRR(find(cumul12on22>0.5,1,'first'));
STATSmodURonRR.cumul = (cumul12on22);
cumul11on21 = cumsum(STATSmodUUonRU.pdf) * (TUUonUR(2)-TUUonUR(1));
STATSmodUUonRU.median = TUUonUR(find(cumul11on21>0.5,1,'first'));
STATSmodUUonRU.mean = ...
    sum(STATSmodUUonRU.pdf .* TUUonUR) * (TUUonUR(2)-TUUonUR(1));
STATSmodUUonRU.cumul = (cumul11on21);
if nargin == 4
    id1         = find(cumul12on22<alphaCI/2,1,'last');
    if isempty(id1)
        STATSmodURonRR.CI(1) = NaN;
    else
        STATSmodURonRR.CI(1) = TURonRR(id1);
    end
    id2         = find(cumul12on22>1-alphaCI/2,1,'first');
    if isempty(id2)
        STATSmodURonRR.CI(2) = NaN;
    else
        STATSmodURonRR.CI(2) = TURonRR(id2);
    end
    
    id1         = find(cumul11on21<alphaCI/2,1,'last');
    if isempty(id1)
        STATSmodUUonRU.CI(1) = NaN;
    else
        STATSmodUUonRU.CI(1) = TUUonUR(id1);
    end
    id2         = find(cumul11on21>1-alphaCI/2,1,'first');
    if isempty(id2)
        STATSmodUUonRU.CI(2) = NaN;
    else
        STATSmodUUonRU.CI(2) = TUUonUR(id2);
    end
end
%=========================================================================
function f = myf(x,xim1,zeta,N)
aux_var = N*log(zeta .* x) -xim1 .* x;
f       = besseli(0,x,1) .* exp(aux_var) ;
%=========================================================================