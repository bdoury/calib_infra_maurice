% calculation of the expectation of Rsup
% as a function of the MSC and the ratio rho of
% the two levels of noises.
%===============================================

clear
addpath  /Users/maurice/etudes/ctbto/allJOBs2015/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/

% global factor for the spectral matrices(no effect)
globalfactor           = 5;

absHUonHR              = 1;
phase_HUonHR_degree    = 18;
argHUonHR_rad          = phase_HUonHR_degree*pi/180;
HUonHR                 = absHUonHR*exp(1j*argHUonHR_rad);
allT.TUUonUR           = linspace(absHUonHR-0.7,absHUonHR+0.7,300);
allT.TURonRR           = linspace(absHUonHR-0.7,absHUonHR+0.7,300);
allT.MSC               = 0;
allT.phase             = linspace(0,360,100);

% C =1/(1+sigmau2)(1+sigmar2)
% sigmar = g x sigmau, with g about 10
% (1+sigmau2)(1+g^2 sigmau2)=1/C
% [1-1/C 1+g^2 g^2]
%
nbnoiseratio = 1;
listnoiseratio = linspace(1
0,10,nbnoiseratio);
nbMSC = 1;
listMSC = linspace(0.96,0.96,nbMSC);

% when it is infinite spectral estimation is
% very accurate

listM = 5:10;
nbM = length(listM);
ERsupmean = zeros(nbM,nbMSC,nbnoiseratio);
ERsupstd = zeros(nbM,nbMSC,nbnoiseratio);


for iM=1:nbM
    twoMmoins1=2*listM(iM)-1;
    for iMSC=1:nbMSC
        
        MSC          = listMSC(iMSC);
        for irho=1:nbnoiseratio
            noiseratio=listnoiseratio(irho);
            sigmau2            = roots([noiseratio^2  ...
                1+noiseratio^2 ...
                1-1/MSC]);
            sigmau2            = max(sigmau2);
            sigmau             = sqrt(sigmau2);
            sigmar             = sigmau*noiseratio;
            sigmar2            = sigmar ^2;
            
            spectralmatrix     = globalfactor * ...
                [(1+sigmar2) / absHUonHR  exp(-1j*argHUonHR_rad); ...
                exp(1j*argHUonHR_rad) (1+sigmau2) * absHUonHR];
            
            [hatpdfUUonUR, hatpdfURonRR, hatMSC] = ...
                statsRatiosHbis(allT,spectralmatrix,twoMmoins1,0.1587);
            ERsupstd(iM,iMSC,irho)=diff(hatpdfUUonUR.CI)/2;
            ERsupmean(iM,iMSC,irho)=(hatpdfUUonUR.mean);
%             plot(allT.TUUonUR,hatpdfUUonUR.cumul)
%             text(1,0.5,num2str(hatpdfUUonUR.mean)) 
%             pause
        end
    end
end

subplot(211)
plot(squeeze(ERsupmean),'.-')




subplot(212)
plot(squeeze(ERsupstd),'.-')






