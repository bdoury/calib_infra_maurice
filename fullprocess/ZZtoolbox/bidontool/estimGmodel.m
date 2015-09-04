function [hatHut, hatgammaSOI, hatsigmau, hatsigmar, fval] = ...
    estimGmodel(Href, xUT, xREF, Lfft, listsigmau, listsigmar)

Lu = length(listsigmau);
Lr = length(listsigmar);
fp = zeros(Lu,Lr);
%==== spectral estimation
Sk11      = mypwelch(xUT,xUT,Lfft,0.5);
Sk22      = mypwelch(xREF,xREF,Lfft,0.5);
Sk12      = mypwelch(xUT,xREF,Lfft,0.5);
for iu=1:Lu,[iu,Lu]
    sigmau2 = listsigmau(iu) ^2;
    for ir =1:Lr,
        sigmar2 = listsigmar(ir) ^2;
        fp(iu,ir) = Jmle_2sigmas(...
            Sk11,...
            Sk22,...
            Sk12,...
            Href,sigmau2,sigmar2);
    end
end

fpmin = min(min(fp));
[indminu,indminr] = find(fp==fpmin);
hatsigmau=listsigmau(indminu);
hatsigmar=listsigmar(indminr);
hatsigmau2 = hatsigmau .^ 2;
hatsigmar2 = hatsigmar .^ 2;

[fp(indminu,indminr),fp(indminr,indminu)]

rhoK = (Sk11 .* Sk22) ./ (abs(Sk12) .^2) - 1;
hatgammaSOI = ((hatsigmau2+hatsigmar2)+ ...
    sqrt((hatsigmau2+hatsigmar2).^2+4 *(hatsigmau2*hatsigmar2) .* rhoK)) ./ ...
    (2*rhoK);
hatHut = Sk12 ./ hatgammaSOI ./ conj(Href);
fval = 0;
for ik=1:Lfft
    Hu = hatHut(ik);
    Hr = Href(ik);
    H = [Hu;Hr];
    HH = H*H';
    Gamma_ik = hatgammaSOI(ik)*HH+diag([hatsigmau2*HH(1,1), hatsigmar2*HH(2,2)]);
    hatGamma = [Sk11(ik) Sk12(ik);conj(Sk12(ik)) Sk22(ik)];
    fval = fval+real(log(det(Gamma_ik)))+real(trace(Gamma_ik\hatGamma));
end
%=======================================================
figure(5)
mesh(listsigmar,listsigmau,fp)
hold on
plot3(listsigmar(indminr),listsigmau(indminu),fpmin,'or','markers',12,'markerfac','r')
title(sprintf('su = %5.2f, sr = %5.2f',hatsigmau,hatsigmar))
hold off

%=======================================================        
function fp = Jmle_2sigmas(Sk11,Sk22,Sk12, Href,sigmau2,sigmar2)


Lfft = length(Href);
gammaK = zeros(Lfft,1);

rhoK = (Sk11 .* Sk22) ./ (abs(Sk12) .^2) - 1;
for ik = 1:Lfft
    rho_ik = rhoK(ik);
    gammaK(ik) = ((sigmau2+sigmar2)+ ...
        sqrt((sigmau2+sigmar2).^2+4 *(sigmau2*sigmar2) .* rho_ik)) ./ ...
        (2*rho_ik);
end
Hut  = (Sk12 ./ gammaK) ./ conj(Href);

fp = 0;
for ik = 1:Lfft
    hatGammaik = [Sk11(ik) Sk12(ik);Sk12(ik)' Sk22(ik)];
    H = [Hut(ik);Href(ik)];
    HHp = H*H';
    Gammatheo = gammaK(ik) * (HHp)+diag([sigmau2*HHp(1,1) sigmar2*HHp(2,2)]);
    fp = fp + real(log(det(Gammatheo)))+real(trace(Gammatheo\hatGammaik));
end

