function [Rsup, Rinf, SCPsupeta, MSC, Nsupthreshold] ...
                 = estimSUT(allSDs, MSCThreshold, flagmedian)
%===========================================================
% synopsis:
% estimate the stastitics of interest from the spectral
%     components. Then extract the values over the 
%     threshold.
%
% [Rsup, Rinf, SpMatrix, MSC, Nsupthreshold] ...
%    = estimSUT(allSDs, MSCThreshold, flagmedian)
%
%====
% Inputs:
%        allSDs:
%            all spectrals components (provided by the function
%                estimSD.
%        MSCThreshold: MSC threshold, typicall 0.99
%        flagmedian: if 1, we take the median, if not we take the mean
%           Rk: the median is not sensible to the outlayers.
%
%====
% Outputs:
%       Rsup: ratio SUU on SUR (see document)
%         structure
%            Rsup.tabmod: table of the modulus of Rsup
%            Rsup.tabmodcst: table of the modulus of Rsup with MSC 
%                   above the threshold
%
%            Rsup.tabphase: table of the phase of Rsup
%            Rsup.tabphasecst: table of the phaseof Rsup with MSC 
%                   above the threshold
%
%            Rsup.mod: modulus of Rsup estimated on Rsup.tabmod
%            Rsup.modcst: modulus of Rsup estimated on Rsup.tabmod
%            Rsup.stdmodcst: STD estimated on Rsup.tabmod 
%                   with MSC above the threshold
5
%            Rsup.phase: table of the phases estimated on Rsup.tabphase
%                   with MSC above the threshold
%            Rsup.phasecst: phases of Rsup estimated on Rsup.tabphase 
%                   with MSC above the threshold;
%            Rsup.stdphasecst: STD estimated on Rsup.tabphase 
%                   with MSC above the threshold
%
%       Rinf: ---- idem Rsup
%
%       SCPsupeta: spectral components 3 x Lfft
%           for cells with MSC is over the threshold. 
%           For each frequency the spectral components is performed
%           as the mean for all cells over the threshold.
%
%       MSC: structure
%           MSC.tab: table of the MSC
%           MSC.tabcst: table of the MSC with values above the threshold
%           MSC.indexcst: index of the 
%               table of the MSC with values above the threshold
%
%       Nsupthreshold: number of values above the threshold
%
%
%===========================================================
%===========================================================

nbblocksAVE = length(allSDs);
Lfft = length(allSDs(1).UU);
indextabMSCsupthreshold = false(Lfft,nbblocksAVE);
indextabMSCsupthreshold([allSDs.MSC]>MSCThreshold) = true;
Nsupthreshold = sum(sum(indextabMSCsupthreshold));

tabUU    = [allSDs.UU];
tabRR    = [allSDs.RR];
tabUR    = [allSDs.UR];
tabMSC   = [allSDs.MSC];

tabMSCwithMSCsupeta = nan(Lfft,nbblocksAVE);
tabMSCwithMSCsupeta(indextabMSCsupthreshold) = ...
    (tabMSC(indextabMSCsupthreshold));

tabUUwithMSCsupeta = nan(Lfft,nbblocksAVE);
tabUUwithMSCsupeta(indextabMSCsupthreshold) = ...
    (tabUU(indextabMSCsupthreshold));

tabRRwithMSCsupeta = nan(Lfft,nbblocksAVE);
tabRRwithMSCsupeta(indextabMSCsupthreshold) = ...
    (tabRR(indextabMSCsupthreshold));
tabURwithMSCsupeta = nan(Lfft,nbblocksAVE);
tabURwithMSCsupeta(indextabMSCsupthreshold) = ...
    (tabUR(indextabMSCsupthreshold));


SCPsupeta      = zeros(3,Lfft);
SCPsupeta(1,:) = nanmean(tabRRwithMSCsupeta,2);
SCPsupeta(2,:) = nanmean(tabUUwithMSCsupeta,2);
SCPsupeta(3,:) = nanmean(tabURwithMSCsupeta,2);

%================================================
%================================================
%================= on Rsup ======================
tabHUUUR    = [allSDs.Rsup];
%=============== to extract the 2D median we take the medians 
%                on x and y axis it is very close to the 
%                minimum of the L1-norm in R2
%===== real/imaginary part without constraint
tabrealHUUUR = real(tabHUUUR);
tabimagHUUUR = imag(tabHUUUR);
tabmodHUUUR  = sqrt(tabrealHUUUR .^2 + tabimagHUUUR .^2);
tabphaseHUUUR = atan2(tabimagHUUUR,tabrealHUUUR);

if flagmedian
        hatrealHUUUR = nanmedian(tabrealHUUUR,2);
        hatimagHUUUR = nanmedian(tabimagHUUUR,2);
else
        hatrealHUUUR = nanmean(tabrealHUUUR,2);
        hatimagHUUUR = nanmean(tabimagHUUUR,2);
end
%===== mod/phase part without constraint
hatmodHUUUR   = sqrt(hatrealHUUUR .^2 + hatimagHUUUR .^2);
hatphaseHUUUR = atan2(hatimagHUUUR, hatrealHUUUR);

%================ on Rsup with constraint ======
%===== real part with constraint
tabrealHUUURwithMSCsupeta = nan(Lfft,nbblocksAVE);
tabrealHUUURwithMSCsupeta(indextabMSCsupthreshold) = ...
    (real(tabHUUUR(indextabMSCsupthreshold)));
%===== imaginary part with constraint
tabimagHUUURwithMSCsupeta = nan(Lfft,nbblocksAVE);
tabimagHUUURwithMSCsupeta(indextabMSCsupthreshold) = ...
    (imag(tabHUUUR(indextabMSCsupthreshold)));
%===== module with constraint
tabmodHUUURwithMSCsupeta = sqrt(...
    tabrealHUUURwithMSCsupeta .^2+...
    tabimagHUUURwithMSCsupeta .^2);
tabphaseHUUURwithMSCsupeta = atan2(tabimagHUUURwithMSCsupeta,tabrealHUUURwithMSCsupeta);

%===== hat's
if flagmedian
    hatrealHUUURwithMSCsupeta = ...
        nanmedian(tabrealHUUURwithMSCsupeta,2);
else
    hatrealHUUURwithMSCsupeta = ...
        nanmean(tabrealHUUURwithMSCsupeta,2);
end
if flagmedian
    hatimagHUUURwithMSCsupeta = ...
        nanmedian(tabimagHUUURwithMSCsupeta,2);
else
    hatimagHUUURwithMSCsupeta = ...
        nanmean(tabimagHUUURwithMSCsupeta,2);
end
%===== mod/phase part with constraint
hatmodHUUURwithMSCsupeta = sqrt(...
    hatrealHUUURwithMSCsupeta .^2+...
    hatimagHUUURwithMSCsupeta .^2);
hatphaseHUUURwithMSCsupeta = atan2(...
    hatimagHUUURwithMSCsupeta,...
    hatrealHUUURwithMSCsupeta);
stdmodHUUURwithMSCsupeta = nanstd(tabmodHUUURwithMSCsupeta,[],2);
stdphaseHUUURwithMSCsupeta = nanstd(tabphaseHUUURwithMSCsupeta,[],2);

%================================================
%================================================
%================= on Rinf ======================
tabHURRR    = [allSDs.Rinf];
%===== real/imaginary part without constraint
tabrealHURRR = real(tabHURRR);
tabimagHURRR = imag(tabHURRR);
tabmodHURRR  = sqrt(tabrealHURRR .^2 + tabimagHURRR .^2);
tabphaseHURRR = atan2(tabimagHURRR, tabrealHURRR);

if flagmedian
        hatrealHURRR = nanmedian(tabrealHURRR,2);
        hatimagHURRR = nanmedian(tabimagHURRR,2);
else
        hatrealHURRR = nanmean(tabrealHURRR,2);
        hatimagHURRR = nanmean(tabimagHURRR,2);
end
%===== mod/phase part without constraint
hatmodHURRR   = sqrt(hatrealHURRR .^2 + hatimagHURRR .^2);
hatphaseHURRR = atan2(hatimagHURRR, hatrealHURRR);

%================ on Rinf with constraint ======
%===== real part with constraint
tabrealHURRRwithMSCsupeta = nan(Lfft,nbblocksAVE);
tabrealHURRRwithMSCsupeta(indextabMSCsupthreshold) = ...
    (real(tabHURRR(indextabMSCsupthreshold)));
%===== imaginary part with constraint
tabimagHURRRwithMSCsupeta = nan(Lfft,nbblocksAVE);
tabimagHURRRwithMSCsupeta(indextabMSCsupthreshold) = ...
    (imag(tabHURRR(indextabMSCsupthreshold)));
%===== module with constraint
tabmodHURRRwithMSCsupeta = sqrt(...
    tabrealHURRRwithMSCsupeta .^2+...
    tabimagHURRRwithMSCsupeta .^2);
tabphaseHURRRwithMSCsupeta = atan2(tabimagHURRRwithMSCsupeta,tabrealHURRRwithMSCsupeta);
%===== hat's
if flagmedian
    hatrealHURRRwithMSCsupeta = ...
        nanmedian(tabrealHURRRwithMSCsupeta,2);
else
    hatrealHURRRwithMSCsupeta = ...
        nanmean(tabrealHURRRwithMSCsupeta,2);
end
if flagmedian
    hatimagHURRRwithMSCsupeta = ...
        nanmedian(tabimagHURRRwithMSCsupeta,2);
else
    hatimagHURRRwithMSCsupeta = ...
        nanmean(tabimagHURRRwithMSCsupeta,2);
end
%===== mod/phase part with constraint
hatmodHURRRwithMSCsupeta = sqrt(...
    hatrealHURRRwithMSCsupeta .^2+...
    hatimagHURRRwithMSCsupeta .^2);
hatphaseHURRRwithMSCsupeta = atan2(...
    hatimagHURRRwithMSCsupeta,...
    hatrealHURRRwithMSCsupeta);

Rsup.tabmod = tabmodHUUUR;
Rsup.tabmodcst = tabmodHUUURwithMSCsupeta;
Rsup.tabphase = tabphaseHUUUR;
Rsup.tabphasecst = tabphaseHUUURwithMSCsupeta;

Rsup.mod = hatmodHUUUR;
Rsup.phase = hatphaseHUUUR;
Rsup.modcst = hatmodHUUURwithMSCsupeta;
Rsup.stdmodcst = stdmodHUUURwithMSCsupeta;
Rsup.phasecst = hatphaseHUUURwithMSCsupeta;
Rsup.stdphasecst = stdphaseHUUURwithMSCsupeta;

Rinf.tabmod = tabmodHURRR;
Rinf.tabmodcst = tabmodHURRRwithMSCsupeta;
Rinf.tabphase = tabphaseHURRR;
Rinf.tabphasecst = tabphaseHURRRwithMSCsupeta;

Rinf.mod = hatmodHURRR;
Rinf.phase = hatphaseHURRR;
Rinf.modcst = hatmodHURRRwithMSCsupeta;
Rinf.phasecst = hatphaseHURRRwithMSCsupeta;

MSC.tab = tabMSC;
MSC.tabcst = tabMSCwithMSCsupeta;
MSC.indexcst = indextabMSCsupthreshold;
%===========================================================
%===========================================================
