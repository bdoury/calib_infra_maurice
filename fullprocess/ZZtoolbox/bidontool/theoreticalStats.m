function [stat11on21, stat12on22, statMSC, STDphase_rd, statPhasepsd] = ...
    theoreticalStats(allT, G, N, alphaCI)

%=========================================================================
% perform
% stat11on21 = statUUonUR, stat12on22 = statURonRR
%=========================================================================
% Compute the probability density functions of
%    - the ratio modules:
%
%             |G(1,2)|                  G(1,1)
%    Rinf =  ----------   and  Rsup =  --------
%              G(2,2)                  |G(2,1)|
%    - the phase
%                phase = arg G(2,1)
%    - the MSC 
%                        |G(1,2)|^2
%                MSC = --------------
%                       G(1,1)G(2,2)
%=========================================================================
% Used Matlab functions
%    INTEGRAL on R2013
%    QUADGK   on R2010
%    besseli
% Inputs :
%       - allT:
%             allT.T11on21: list of values range of S11/abs(S21)
%             allT.T12on22: list of values range of S12/abs(S22)
%             allT.phase: list of values of arg(S22)
%             allT.MSC: list of values of MSC
%       - G: is the 2 x 2 spectral matrix at one frequency bin
%                   | G(1,1) G(1,2) |
%                   | G(2,1) G(2,2) |
%             G(1,1) is spectrum on SUT, G(2,2) is spectrum on SREF
%       - N: window length for averaging FFT frames
%       - alphaCI: level of confidence, typical values 0.05, 0.15, 0.3
%=========================================================================
allT.T11on21 = allT.TUUonUR;
allT.T12on22 = allT.TURonRR;
% we need to re-arrange G
G = [G(2,2) G(2,1);G(1,2) G(1,1)];
% values of MSC range
valMSC              = allT.MSC;
%=======
rho                 = abs(G(1,2)/sqrt(G(1,1)*G(2,2)));
rho2                = rho ^2;
lambda              = sqrt(G(2,2)/G(1,1));
% (1*2*...*(N-1))^(1/N)
gammaNm1P           = exp(sum(log(1:N-1))/N);
MSC_theo            = abs(G(2,1)) ^2/G(1,1)/G(2,2);
phase_HUminusHR_rad = -atan2(imag(G(2,1)),real(G(2,1)));
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
T         = allT.T12on22;
T         = T(:);
zeta      = (1-rho2) ./ (2*rho*lambda * T)/gammaNm1P;
cst       = lambda/rho;
xi        = ((1+lambda*lambda*(T .* T))) ./(2*rho*lambda*T);
xim1      = xi-1;
LT        = length(T);
p12on22   = zeros(LT,1);
for it=1:LT
    if exist('integral','file')
        p12on22(it)   = cst * integral(@(x) myf(x,xim1(it), ...
            zeta(it),N),0,inf);
    else
        p12on22(it)   = cst * quadgk(@(x) myf(x,xim1(it), ...
            zeta(it),N),0,inf);
    end
end
stat12on22.pdf = p12on22;
valMSC         = valMSC(:);
statMSC.pdf    = pdfMSC(valMSC,rho2,N);
statMSC.CI(1)  = invcumulFunctionMSC(alphaCI/2,rho2,N);
statMSC.CI(2)  = invcumulFunctionMSC(1-alphaCI/2,rho2,N);
%======================================================================
%===================== second ratio ===================================
%======================================================================
Tper      = 1 ./ allT.T11on21;
Tper      = Tper(:);
Gper      = [G(2,2) G(2,1);G(1,2) G(1,1)];
rho       = abs(Gper(1,2)/sqrt(Gper(1,1)*Gper(2,2)));
rho2      = rho ^2;
lambda    = sqrt(Gper(2,2)/Gper(1,1));
zeta      = (1-rho2) ./ (2*rho*lambda * Tper)/gammaNm1P;
cst       = lambda/rho;
xi        = ((1+lambda*lambda*(Tper .* Tper)))...
    ./(2*rho*lambda*Tper);
xim1      = xi-1;
LTper     = length(Tper);
p11on21   = zeros(LTper,1);
for it=1:LTper
    if exist('integral','file')
        p11on21(it)   = cst * integral(@(x) myf(x,xim1(it), ...
            zeta(it),N),0,inf);
    else
        p11on21(it)   = cst * quadgk(@(x) myf(x,xim1(it), ...
            zeta(it),N),0,inf);
    end
end
stat11on21.pdf = p11on21 .* (Tper .* Tper) ;
cumul12on22 = cumsum(stat12on22.pdf)*...
    (allT.T12on22(2)-allT.T12on22(1));
stat12on22.mean = ...
    sum(stat12on22.pdf .* allT.T12on22')*...
    (allT.T12on22(2)-allT.T12on22(1));
stat12on22.median = allT.T12on22(find(cumul12on22>0.5,1,'first'));
stat12on22.cumul = (cumul12on22);
cumul11on21 = cumsum(stat11on21.pdf)*...
    (allT.T11on21(2)-allT.T11on21(1));
stat11on21.median = allT.T11on21(find(cumul11on21>0.5,1,'first'));
stat11on21.mean = ...
    sum(stat11on21.pdf .* allT.T11on21')*...
    (allT.T11on21(2)-allT.T11on21(1));
stat11on21.cumul = (cumul11on21);
if nargin == 4
    id1         = find(cumul12on22<alphaCI/2,1,'last');
    if isempty(id1)
        stat12on22.CI(1) = NaN;
    else
        stat12on22.CI(1) = allT.T12on22(id1);
    end
    id2         = find(cumul12on22>1-alphaCI/2,1,'first');
    if isempty(id2)
        stat12on22.CI(2) = NaN;
    else
        stat12on22.CI(2) = allT.T12on22(id2);
    end
    
    id1         = find(cumul11on21<alphaCI/2,1,'last');
    if isempty(id1)
        stat11on21.CI(1) = NaN;
    else
        stat11on21.CI(1) = allT.T11on21(id1);
    end
    id2         = find(cumul11on21>1-alphaCI/2,1,'first');
    if isempty(id2)
        stat11on21.CI(2) = NaN;
    else
        stat11on21.CI(2) = allT.T11on21(id2);
    end
end
%=========================================================================
function f = myf(x,xim1,zeta,N)
aux_var = N*log(zeta .* x) -xim1 .* x;
f       = besseli(0,x,1) .* exp(aux_var) ;
%=========================================================================