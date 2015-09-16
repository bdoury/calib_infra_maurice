function  [SUTs, filteredsignals, allfrqsFFT_Hz, alltimes_sec, filterbank] = ...
    fbankanalysis(sigin, ...
    filtercharacteristics, ...
    Fs_Hz,...
    MSCthreshold)
%===========================================================
% synopsis:
%    SUTs = fbankanalysis...
%          (sigin,filtercharacteristics,Fs_Hz,MSCthreshold)
% DFT for discrete Fourier transform
% SCP for spectral components
%     SCP_k = MEAN_{m=1}^M wDFT_{m,k}
%     wDFT_{m,k} = DFT( window .* signal_m, k, Lfft)
%          where signal_m is the m-th segment with overlap,
%          Lfft the DFT length
%          and k the frequency index associated to the
%          frequency k*Fs_Hz/Lfft
%===========================================================
% Remark:
% The 3 parameters, M, Lfft, and stationarity-duration are
% related by:
%              stationarity-duration ~ M * Lfft / Fs_Hz
%
% The choice of M is related to the accuracy of the estimator
% In the following we focus on the stationarity-duration value,
% depending on the frequency band, and on M:
%    stationarity-duration is defined by xx.SCPperiod_sec
%    M is defined by xx.ratioDFT2SCP (overlap not applied)
%    therefore we derive Lfft, regarding the overlap
%    xx.overlapDFT
%===========================================================
% Inputs:
%    sigin: input signals T x 2
%    filtercharacteristics: structure Px1
%         P: number of frequency bands
%         xx.Norder: order of the filter
%                if xx.Norder=0, no filtering
%         xx.Wlow_Hz: lower bound of the frequency band in Hz
%         xx.Whigh_Hz: upper bound of the frequency band in Hz
%         xx.SCPperiod_sec: duration in second
%                 over which SPC is performed,
%                 expected as stationarity time
%         xx.windowshape: window shape for spectral
%                 analysis (ex. 'hann' or 'hamming', ...)
%         xx.overlapDFT: overlap rate for successive
%                 DFT (typically 0.5)
%         xx.overlapSCP: overlap rate for successive
%                 spectral components (typically 0)
%         xx.ratioDFT2SCP: ratio between period_sec and
%                 DFT duration (typical integer value is 5).
%     Fs_Hz: sampling frequency in Hz
%     MSCthreshold: MSC threshold (advised value > 0.99)
%
%===========================================================
% Rk: the spectral components are performed on successive DFTs with
%     overlapping. If M denotes the ratioFFT2DSP and
%     if overlapFFT is 0.5, the number of DFTs is 9
%===========================================================
% Example of filtercharacteristics structure:
%===========================================================
%         xx.Norder: 4
%         xx.Wlow_Hz: 0.02
%         xx.Whigh_Hz: 0.2
%         xx.SCPperiod_sec: 250
%         xx.windowshape: 'hann'
%         xx.overlapDFT: 0.5
%         xx.overlapSCP: 0
%         xx.ratioDFT2SCP: 5.
%===========================================================
% filtercharact(1).Norder           = 4;
% filtercharact(1).Wlow_Hz          = 0.02;
% filtercharact(1).Whigh_Hz         = 0.2;
% filtercharact(1).SCPperiod_sec    = 500;
% filtercharact(1).windowshape      = 'hann';
% filtercharact(1).overlapDFT       = 0.5;
% filtercharact(1).overlapSCP       = 0;
% filtercharact(1).ratioDFT2SCP     = 5;
%===========================================================
%===========================================================
% Outputs:
%     SUTs: structures P x 1
%         xx.estimRsup: Rsup ratio
%         xx.estimRinf: Rinf ratio
%         xx.allMSCs: allMSCs;
%         xx.Nsupthreshold: count of the number of values over
%                   the threshold
%         xx.Nsupthresholdintheband: counts of the number of values
%                   over the threshold in the filter bandwidth
%         xx.frqsFFT_Hz: cell P x 1, each cell consists of
%                   frequency list in Hz of the DFTs
%         xx.SCP = all spectral components
%         xx.indexinsidefreqband = P x 1, indices of the
%                    frequency bounds of each filter
%                    in the xx.frqsFFT_Hz
%         xx.alltimes_sec: cell  Px 1, each cell consists of
%             yy.FFT: time list in second of the DFTs
%             yy.SD: time list in second of the SCPs
%             yy.signals: time list in second of the signals
%===========================================================
% included functions: estimSCPs, estimSUT, fbank
%===========================================================
%===========================================================
overlapDFT       = filtercharacteristics.overlapDFT;
overlapSCP       = filtercharacteristics.overlapSCP;
ratioDFT2SCP     = filtercharacteristics.ratioDFT2SCP;
P                = length(filtercharacteristics);
allfrqsFFT_Hz    = cell(P,1);
alltimes_sec     = cell(P,1);
SUTs             = struct;
T                = size(sigin,1);
%===== filtering SIGIN ==> SIGOUT
[sigout, filterbank]  = fbank(sigin,filtercharacteristics,Fs_Hz);
filteredsignals       = zeros(T,2,P);
for ifilter = 1:P
    taustationary_sec = filtercharacteristics(ifilter).SCPperiod_sec;
    windowingshape    = filtercharacteristics(ifilter).windowshape;
    signals   = sigout(:,:,ifilter);
    [T,d]     = size(signals);
    signals   = signals-ones(T,1)*mean(signals);
    filteredsignals(:,:,ifilter) = signals;
    Tfft_sec  = taustationary_sec/ratioDFT2SCP;
    Lfft      = fix(Tfft_sec*Fs_Hz);
    %===================== spectral analysis =========================
    [allSDs, time_temp, allfrqsFFT_Hz{ifilter}] = ...
        estimSCP(signals(:,1),signals(:,2),...
        ones(Lfft,1), overlapDFT, ...
        ratioDFT2SCP,overlapSCP, Fs_Hz, windowingshape);
    alltimes_sec{ifilter} = time_temp;
    %================================================================
    [estimRsup, estimRinf, spectralmatrix, allMSCs, Nsupthreshold] ...
        = estimSUT(allSDs, MSCthreshold);
    SUTs(ifilter).estimRsup = estimRsup;
    SUTs(ifilter).estimRinf = estimRinf;
    SUTs(ifilter).allMSCs   = allMSCs;
    SUTs(ifilter).Nsupthreshold = Nsupthreshold;
    SUTs(ifilter).frqsFFT_Hz = allfrqsFFT_Hz{ifilter};
    SUTs(ifilter).SCP = allSDs;
    SUTs(ifilter).indexinsidefreqband = [ ...
        ceil((filtercharacteristics(ifilter).Wlow_Hz/Fs_Hz)*length(SUTs(ifilter).frqsFFT_Hz)); ...
        ceil((filtercharacteristics(ifilter).Whigh_Hz/Fs_Hz)*length(SUTs(ifilter).frqsFFT_Hz))];
    if and(P==1, filtercharacteristics(ifilter).Norder==0)
                SUTs(ifilter).Nsupthresholdintheband = ...
            nansum(SUTs(ifilter).allMSCs.indexcst,2);
    else
        SUTs(ifilter).Nsupthresholdintheband = ...
            nansum(SUTs(ifilter).allMSCs.indexcst(SUTs(ifilter).indexinsidefreqband(1):...
            SUTs(ifilter).indexinsidefreqband(2),:),2);
    end
end
%=======================================================
function [sigout, filterbank] = fbank(sigin,filtercharact,Fs_Hz)
%=======================================================
P           = length(filtercharact);
[T,d]       = size(sigin);
sigout      = zeros(T,d,P);
filterbank = cell(P,1);
for ifilter = 1:P
    if not(filtercharact(ifilter).Norder==0)
        flow = filtercharact(ifilter).Wlow_Hz/Fs_Hz;
        fhigh = filtercharact(ifilter).Whigh_Hz/Fs_Hz;
        fname = filtercharact(ifilter).designname;
        forder = filtercharact(ifilter).Norder;
        switch fname
            case 'fir1'
                fdesign = sprintf('filnum = %s(%i,[%5.8f,%5.8f]);',...
                    fname,forder,2*flow,2*fhigh);
                filden = 1;
            case 'butter'
                fdesign = sprintf('[filnum,filden] = %s(%i,[%5.8f %5.8f]);',...
                    fname,forder,2*flow,2*fhigh);
            case 'cheby1'
                fdesign = sprintf('[filnum,filden] = %s(%i,%i,[%5.8f %5.8f]);',...
                    fname,forder,0.02,2*flow,2*fhigh);
        end
        eval(fdesign);
        sigout(:,:,ifilter) = filter(filnum,filden,sigin);
        filterbank{ifilter}.num = filnum;
        filterbank{ifilter}.den = filden;
        
    else
        sigout(:,:,ifilter) = sigin;
        filterbank{ifilter}.num = NaN;
        filterbank{ifilter}.den = NaN;
        
    end
end
%=======================================================
function [allSDs, time_sec, frqsFFT_Hz] = ...
    estimSCP(xU,xR,GREF,overlapFFT, ...
    NaverageFFTs, overlapSD, Fs_Hz, smoothwindow)
%==========================================================================
% Perform the spectral components of the two signals xU et xR.
%==========================================================================
% The code uses the Welch's approach. The signal is shared into DFT
% windows of which the length is the length of GREF, with the
% overlap rate of OVERLAPFFT. Then the specral components is averaged
% on NaverageFFTs DFT blocks. Therefore each spectral block corresponds
% to a time period reported in TIME_SEC.
%
% Inputs:
%    xU: signal observed on the SUT (T x 1)
%    xR: signal observed on the SREF (T x 1)
%    GREF: frequency response of the SREF (Lfft x 1)
%    overlapFFT: between 0 and 1, overlap on the FFT block
%    NaverageFFTs: number of FFT-length to averaging
%    overlapSD: between 0 and 1, overlap on the averaging
%                  Spectral Density block
%    Fs_Hz: sampling frequency in Hz
%    smoothwindow: character word as 'rect', 'hann', ...
%            [default: rectangular]
% Outputs:
%    allSDs.RR = auto-spectrum of SREF
%    allSDs.UU = auto-spectrum of SUT
%    allSDs.UR = cross-spectrum of (SUT,SREF)
%    allSDs.MSC = Magnitude Square Coherence
%    allSDs.Rsup = SUU/SUR
%    allSDs.Rinf = SUR/SRR
%    allSDs.det = determinant of the spectral matrix
%
%    time_sec.FFT: grid of time for the FFT blocks (in second)
%    time_sec.SD: grid of time for the Spectral Density blocks (in second)
%    time_sec.signals: grid of time for the samplig period (in second)
%    frqsFFT_Hz: grid of frequency for the interval [0, Fs_Hz] (in Hz)
%==========================================================================
%==========================================================================
xU   = xU(:);
xR   = xR(:);
N    = length(xU);
Lfft = length(GREF);
sqrtLfft    = sqrt(Lfft);
shiftSignal = fix((1-overlapFFT)*Lfft);
NblocksFFT  = fix((N-(Lfft-shiftSignal))/shiftSignal);
allFFTsRR   = zeros(Lfft,NblocksFFT);
allFFTsUU   = zeros(Lfft,NblocksFFT);
if exist('smoothwindow','var')
    switch smoothwindow
        case 'hann'
            Hwin = hann(Lfft,'periodic');
        case 'bartlett'
            Hwin = bartlett(Lfft);
        case 'hamming'
            Hwin = hamming(Lfft);
        case 'rectwin'
            Hwin = rectwin(Lfft);
        case 'blackman'
            Hwin = blackman(Lfft);
    end
else
    Hwin     = hamming(Lfft);
end
%========= normalisation
% not useful if only PSD ratios are considered
Hwin = Hwin *sqrt(Lfft/(Hwin'*Hwin));
for ibF  = 1:NblocksFFT
    ibT  = (ibF-1)*shiftSignal+(1:Lfft);
    xU_i = xU(ibT) .* Hwin;
    xU_i = xU_i-mean(xU_i);
    xR_i = xR(ibT) .* Hwin;
    xR_i = xR_i-mean(xR_i);
    allFFTsUU(:,ibF) = fft(xU_i,Lfft)/sqrtLfft;
    allFFTsRR(:,ibF) = fft(xR_i,Lfft)/sqrtLfft;
end
shiftFFTs = fix((1-overlapSD)*NaverageFFTs);
NSD       = fix(N/Lfft/NaverageFFTs);% fix((NblocksFFT/2-(NaverageFFTs-shiftFFTs))/shiftFFTs);
allSDs    = struct;
time_sec  = struct;
NaverageFFTs_with_overlap = round(NaverageFFTs/(1-overlapFFT)-1);
for ibB=1:NSD,
    indB1 = (ibB-1)*shiftFFTs+1;
    indB2 = indB1+NaverageFFTs_with_overlap-1;
    indB  = fix(indB1):fix(indB2);
    indB  = indB(indB<= NblocksFFT);
    allSDs(ibB).RR  = mean(abs(allFFTsRR(:,indB)) .^ 2,2);
    allSDs(ibB).UU  = mean(abs(allFFTsUU(:,indB)) .^ 2,2);
    allSDs(ibB).UR  = mean(allFFTsRR(:,indB) .* conj(allFFTsUU(:,indB)),2);
    allSDs(ibB).MSC = (abs(allSDs(ibB).UR) .^2) ./...
        (allSDs(ibB).RR .* allSDs(ibB).UU);
    allSDs(ibB).Rsup  = allSDs(ibB).UU ./ allSDs(ibB).UR;
    allSDs(ibB).Rinf  = allSDs(ibB).UR ./ allSDs(ibB).RR;
    allSDs(ibB).det   = (abs(allSDs(ibB).RR) .* abs(allSDs(ibB).UU)) ...
        -(abs(allSDs(ibB).UR) .^2);
end
frqsFFT_Hz = (0:Lfft-1)*Fs_Hz/Lfft;
time_sec.FFT = ((0:NblocksFFT-1)+1/2)*shiftSignal/Fs_Hz;
time_sec.SD =  ((0:NSD-1)+1/2)*shiftFFTs*Lfft/Fs_Hz;
time_sec.signals =  ((0:N-1)+1/2)/Fs_Hz;
% %==========================================================================
function [Rsup, Rinf, SCPsupeta, MSC, Nsupthreshold] ...
    = estimSUT(allSDs, MSCThreshold)
%===========================================================
% synopsis:
% estimate the stastitics of interest from the spectral
%     components. Then extract the values over the
%     threshold.
%
% [Rsup, Rinf, SpMatrix, MSC, Nsupthreshold] ...
%    = estimSUT(allSDs, MSCThreshold)
%
%====
% Inputs:
%        allSDs:
%            all spectrals components (provided by the function
%                estimSD.
%        MSCThreshold: MSC threshold, typicall 0.99
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

%================================================
%================ if weightingflag = 1 ==========
% a ponderation is applied using the estimated variance
% of the MSC estimate
weightingflag = 1;

nbblocksAVE   = length(allSDs);
Lfft          = length(allSDs(1).UU);
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
%================= on Rsup, Rinf======================
tabHUUUR    = [allSDs.Rsup];
tabHURRR    = [allSDs.Rinf];
%================================================
%================================================
%================= on Rsup ======================
%===== real/imaginary part without constraint
tabrealHUUUR = real(tabHUUUR);
tabimagHUUUR = imag(tabHUUUR);
tabmodHUUUR  = sqrt(tabrealHUUUR .^2 + tabimagHUUUR .^2);
tabphaseHUUUR = atan2(tabimagHUUUR,tabrealHUUUR);

%================================================
%================================================
%================= on Rinf ======================
%===== real/imaginary part without constraint
tabrealHURRR = real(tabHURRR);
tabimagHURRR = imag(tabHURRR);
tabmodHURRR  = sqrt(tabrealHURRR .^2 + tabimagHURRR .^2);
tabphaseHURRR = atan2(tabimagHURRR, tabrealHURRR);

%================================================
%================================================
%================= on Rsup HAT ==================
hatrealHUUUR = nanmean(tabrealHUUUR,2);
hatimagHUUUR = nanmean(tabimagHUUUR,2);
%===== mod/phase part without constraint
hatmodHUUUR   = sqrt(hatrealHUUUR .^2 + hatimagHUUUR .^2);
hatphaseHUUUR = atan2(hatimagHUUUR, hatrealHUUUR);
%================= on Rinf HAT ==================
hatrealHURRR = nanmean(tabrealHURRR,2);
hatimagHURRR = nanmean(tabimagHURRR,2);
%===== mod/phase part without constraint
hatmodHURRR   = sqrt(hatrealHURRR .^2 + hatimagHURRR .^2);
hatphaseHURRR = atan2(hatimagHURRR, hatrealHURRR);

%================================================
%================================================
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


%================================================
%================================================
%================ on Rsinf with constraint ======
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

RsuptheowithMSCsupeta = tabmodHUUURwithMSCsupeta .* ...
    tabmodHURRRwithMSCsupeta;

%================================================
%================ if weightingflag = 1 ==========
% a ponderation is applied using the estimated variance
% of the MSC estimate
weightMSCsupeta = (tabMSCwithMSCsupeta .^2) ./  ...
    (1-tabMSCwithMSCsupeta) .* RsuptheowithMSCsupeta;
weightMSCinfeta = 1 ./  ...
    (1-tabMSCwithMSCsupeta) .* RsuptheowithMSCsupeta;

%================================================
%================================================
%========== on Rsup HAT with constraint =========
%========== and weighting coeffs ================

if weightingflag
    hatrealHUUURwithMSCsupeta = ...
        nansum(tabrealHUUURwithMSCsupeta .* weightMSCsupeta,2) ...
        ./ nansum((weightMSCsupeta),2);
    
    hatimagHUUURwithMSCsupeta = ...
        nansum(tabimagHUUURwithMSCsupeta .* weightMSCsupeta,2) ...
        ./ nansum((weightMSCsupeta),2);
else
    hatrealHUUURwithMSCsupeta = ...
        nanmean(tabrealHUUURwithMSCsupeta,2);
    
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

stdmodHUUURwithMSCsupeta   = nanstd(tabmodHUUURwithMSCsupeta,[],2);
stdphaseHUUURwithMSCsupeta = nanstd(tabphaseHUUURwithMSCsupeta,[],2);

%================================================
%================================================
%========== on Rinf HAT with constraint =========
%========== and weighting coeffs ================

if weightingflag
    hatrealHURRRwithMSCsupeta = ...
        nansum(tabrealHURRRwithMSCsupeta .* weightMSCinfeta,2) ...
        ./ nansum((weightMSCinfeta),2);
    
    hatimagHURRRwithMSCsupeta = ...
        nansum(tabimagHURRRwithMSCsupeta .* weightMSCinfeta,2) ...
        ./ nansum((weightMSCinfeta),2);
else
    hatrealHURRRwithMSCsupeta = ...
        nanmean(tabrealHURRRwithMSCsupeta,2);
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

stdmodHURRRwithMSCsupeta   = nanstd(tabmodHURRRwithMSCsupeta,[],2);
stdphaseHURRRwithMSCsupeta = nanstd(tabphaseHURRRwithMSCsupeta,[],2);


%================ on Rinf with constraint ======
%===== hat's

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
Rinf.stdmodcst = stdmodHURRRwithMSCsupeta;
Rinf.phasecst = hatphaseHURRRwithMSCsupeta;
Rinf.stdphasecst = stdphaseHURRRwithMSCsupeta;

MSC.tab = tabMSC;
MSC.tabcst = tabMSCwithMSCsupeta;
MSC.indexcst = indextabMSCsupthreshold;
MSC.weightMSC = weightMSCsupeta;
%===========================================================
%===========================================================


