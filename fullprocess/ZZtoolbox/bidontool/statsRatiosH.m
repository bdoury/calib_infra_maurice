function [stat12on22, stat11on21, statMSC] = ...
    statsRatiosH(allT, G, N, a)
%===============================================
% Compute the probabiity density function of
% the ratios:
%
%       |G(1,2)|             G(1,1)
%      ----------   and    ----------
%        G(2,2)             |G(2,1)|
% where G is a positive matrix.
%
% Derive the confidence interval at 100alpha%
% To  integrate use:
%    INTEGRAL on R2013
%    QUADGK   on R2010
% 
% Inputs : 
%       - allT:
%       - G: spectral matrix at one frequency bin
%       - N: window length for averaging FFT frames
%       - a: level of confidence, typically 0.05%
%======================================================================
% 1=U, R=2
allT.T12on22 = allT.TURonRR;
allT.T11on21 = allT.TUUonUR; 

rho = abs(G(1,2)/sqrt(G(1,1)*G(2,2)));
rho2 = rho ^2;
lambda = sqrt(G(2,2)/G(1,1));
% (1*2*...*(N-1))^(1/N)
gammaNm1P = exp(sum(log(1:N-1))/N);
%======================================================================
%====================== first ratio ===================================
%======================================================================
T = allT.T12on22;
T = T(:);

zeta = (1-rho2) ./ (2*rho*lambda * T)/gammaNm1P;
cst  = lambda/rho;

xi = ((1+lambda*lambda*(T .* T))) ./(2*rho*lambda*T);
xim1 = xi-1;

LT = length(T);
p12on22 = zeros(LT,1);

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


E = allT.MSC;
E = E(:);
statMSC.pdf   = pdfMSC(E,rho2,N);
statMSC.CI(1) = invcumulFunctionMSC(a/2,rho2,N);
statMSC.CI(2) = invcumulFunctionMSC(1-a/2,rho2,N);

%======================================================================
%===================== second ratio ===================================
%======================================================================
Tper = 1 ./ allT.T11on21;
Tper = Tper(:);
Gper = [G(2,2) G(2,1);G(1,2) G(1,1)];
rho = abs(Gper(1,2)/sqrt(Gper(1,1)*Gper(2,2)));
rho2 = rho ^2;
lambda = sqrt(Gper(2,2)/Gper(1,1));


zeta = (1-rho2) ./ (2*rho*lambda * Tper)/gammaNm1P;
cst  = lambda/rho;

xi = ((1+lambda*lambda*(Tper .* Tper))) ./(2*rho*lambda*Tper);
xim1 = xi-1;

LTper = length(Tper);
p11on21 = zeros(LTper,1);

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
    
    id1         = find(cumul12on22<a/2,1,'last');
    if isempty(id1)
        stat12on22.CI(1) = NaN;
    else
        stat12on22.CI(1) = allT.T12on22(id1);
    end
    id2         = find(cumul12on22>1-a/2,1,'first');
    if isempty(id2)
        stat12on22.CI(2) = NaN;
    else
        stat12on22.CI(2) = allT.T12on22(id2);
    end
    
    id1         = find(cumul11on21<a/2,1,'last');
    if isempty(id1)
        stat11on21.CI(1) = NaN;
    else
        stat11on21.CI(1) = allT.T11on21(id1);
    end
    id2         = find(cumul11on21>1-a/2,1,'first');
    if isempty(id2)
        stat11on21.CI(2) = NaN;
    else
        stat11on21.CI(2) = allT.T11on21(id2);
    end
end
%===============================================
function f = myf(x,xim1,zeta,N)
aux_var = N*log(zeta .* x) -xim1 .* x;
f = besseli(0,x,1) .* exp(aux_var) ;
%===============================================



