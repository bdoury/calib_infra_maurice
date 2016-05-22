function theta = optimWhittle(Sigma, Gref, thetainit)

[d,d,K] = size(Sigma);
dKp2 = length(thetainit);
theta = fmincon(@(x) Jof(x,Sigma,Gref), thetainit,[],[],[],[],0.1*ones(dKp2,1));


function J=Jof(x,Sigma,Gref)
[d,d,K] = size(Sigma);
gamma = x(1:K);
Gu = x(K+1:2*K);
s1=x(2*K+1);
s2=x(2*K+2);
J=0;
for ik = 1:K,
    Rhat = squeeze(Sigma(:,:,ik));
    Rik=[Gref(ik);Gu(ik)]*[Gref(ik);Gu(ik)]'*gamma(ik)+diag([s1;s2]);
    J = J+real(log(det(Rik))+trace(Rik\Rhat));
end