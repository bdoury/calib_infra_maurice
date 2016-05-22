%====================================================
% RMSE by calculation on synthetic data
%====================================================
% Perform for a given value of the ratio of HUT on HREF 
% and a given value of the size of the averaging window
% the RMSE.
% This program is called by MULTIRUNS which saves the
% values on .mat file
% Then the program ALLPRINTSTATONHEST.m is called.
%====================================================
% In this case the 2 following rows must be comment.
%====================================================
clear 
%====================================================
addpath  /Users/maurice/etudes/ctbto/allJOBs2015/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/
%====================================================
% R=1, U=2
% we loock for statistics of ratioUUonUR, ratioURonRR, 
% then HU = HR x ratioUUonUR OR HU = HR x ratioURonRR
%====================================================
absHUonHR              = 1;
% global factor for the spectral matrices(no effect)
globalfactor           = 5;
phase_HUonHR_degree    = 18;
argHUonHR_rad          = phase_HUonHR_degree*pi/180;
HUonHR                 = absHUonHR*exp(1j*argHUonHR_rad);
allT.TUUonUR           = linspace(absHUonHR-0.7,absHUonHR+0.7,300);
allT.TURonRR           = linspace(absHUonHR-0.6,absHUonHR+0.7,300);
allT.MSC               = 0;
allT.phase             = linspace(0,360,100);

% C =1/(1+sigmau2)(1+sigmar2)
% sigmar = g x sigmau, with g about 10
% (1+sigmau2)(1+g^2 sigmau2)=1/C
% [1-1/C 1+g^2 g^2]
%
holenumber = 68;
noiseratio = sqrt(holenumber);

listtwoMmoins1         = 10;
L2Mm1                  = length(listtwoMmoins1);
listMSC_theo           = 0.8:0.02:0.99;
LMSC                   = length(listMSC_theo);
stdonHUonHRex1         = zeros(LMSC,L2Mm1);
stdonHUonHRex2         = zeros(LMSC,L2Mm1);

for it = 1:L2Mm1,it
    twoMplus1              = listtwoMmoins1(it);
    for ir = 1:LMSC
        MSC_theo           = listMSC_theo(ir);
        sigmau2            = roots([noiseratio^2  ...
            1+noiseratio^2 ...
            1-1/MSC_theo]);
        sigmau2            = max(sigmau2);
        sigmau             = sqrt(sigmau2);
        sigmar             = sigmau*noiseratio;
        sigmar2            = sigmar ^2;

        spectralmatrix     = globalfactor * ...
            [(1+sigmar2) / absHUonHR  exp(-1j*argHUonHR_rad); ...
            exp(1j*argHUonHR_rad) (1+sigmau2) * absHUonHR];

        [hatpdfUUonUR, hatpdfURonRR, hatMSC] = ...
            statsRatiosHbis(allT,spectralmatrix,twoMplus1,0.05);
        
        stdonHUonHRex1(ir,it) = ...
            sqrt((absHUonHR-hatpdfUUonUR.mean)^2+(diff(hatpdfUUonUR.CI)/2)^2);
        
        stdonHUonHRex2(ir,it) = ...
            sqrt((absHUonHR-hatpdfURonRR.mean)^2+(diff(hatpdfURonRR.CI)/2)^2);
    end
end
save statsonHest1.mat
%=======================================
