function [NWstat,testStat, H] = kpsstat(y, Tc)

alphagrid = [0.010 0.025 0.050 0.100];
Cofalpha  = [0.216  0.176  0.146  0.119];
% Cofalpha  = [0.739  0.574  0.463  0.347];
confidenceLevel = 0.05;

y=y(:);
y(isnan(y))=[];
T = length(y);
H = [ones(T,1) (1:T)'];
beta = H\y;
yhat = H*beta;
res = y - yhat;

gamma = zeros(Tc+1,1);
for k=1:Tc+1
    gamma(k)=res(1:T-k+1)'*res(k:T)/T;
end
S=0;
for k=1:Tc
    S = S + (1-k/(1+Tc))*gamma(k+1);
end
NWstat = gamma(1)+2*S;

eSum = cumsum(res);
testStat = (eSum'*eSum)/(NWstat*T^2);
interpinvcumul = interp1(alphagrid,Cofalpha,confidenceLevel,'linear');
H = (testStat > interpinvcumul); 