%========================================================================
function [SUTs, filteredsignals, allfrqsFFT_Hz, alltimes_sec, ...
    filterbank] = fbankanalysis(...
    sigin, ...
    filtercharacteristics, ...
    Fs_Hz,...
    MSCthreshold,...
    flagtheoreticalSTDs)
%========================================================================
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
%========================================================================
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
%========================================================================
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
%                 DFT duration (typical integer value is 5)
%     Fs_Hz: sampling frequency in Hz
%     MSCthreshold: MSC threshold (advised value > 0.98)
%     flagtheoreticalSTDs: if 1, the theoretical STDs 
%         (fields xx.theoSTDmodonRsup and theoSTDphase_rd of SUT variable) 
%         are performed
%========================================================================
% Rk: the spectral components are performed on successive DFTs with
%     overlapping. If M denotes the ratioFFT2DSP and
%     if overlapFFT is 0.5, the number of DFTs is 9
%===========================================================
% Example of filtercharacteristics structure:
%========================================================================
%         xx.Norder: 4
%         xx.Wlow_Hz: 0.02
%         xx.Whigh_Hz: 0.2
%         xx.SCPperiod_sec: 250
%         xx.windowshape: 'hann'
%         xx.overlapDFT: 0.5
%         xx.overlapSCP: 0
%         xx.ratioDFT2SCP: 5.
%========================================================================
% filtercharact(1).Norder           = 4;
% filtercharact(1).Wlow_Hz          = 0.02;
% filtercharact(1).Whigh_Hz         = 0.2;
% filtercharact(1).SCPperiod_sec    = 500;
% filtercharact(1).windowshape      = 'hann';
% filtercharact(1).overlapDFT       = 0.5;
% filtercharact(1).overlapSCP       = 0;
% filtercharact(1).ratioDFT2SCP     = 5;
%========================================================================
%========================================================================
% Outputs:
%     SUTs: structures P x 1
%         xx.estimRsup: Rsup ratio
%         xx.estimRinf: Rinf ratio
%         xx.allMSCs: allMSCs;
%         xx.Nsupthreshold: count of the number of values over
%                   the threshold
%         xx.Nsupthresholdintheband: counts of the number of values
%                   over the threshold in the filter bandwidth
%         xx.SCP = all spectral components
%         xx.indexinsidefreqband = P x 1, indices of the
%                    frequency bounds of each filter
%                    in the xx.frqsFFT_Hz
%         xx.theoSTDmodonRsup: theoretical STD of the module estimate
%         xx.theoSTDphase_rd: theoretical STD of the phase estimate
%
%         alltimes_sec: cell  Px 1, each cell consists of
%             yy.FFT: time list in second of the DFTs
%             yy.SD: time list in second of the SCPs
%             yy.signals: time list in second of the signals
%
%         frqsFFT_Hz: cell P x 1, each cell consists of
%                   frequency list in Hz of the DFTs
%========================================================================
% included functions: estimSCP, estimSUT, fbank, theoreticalStats
%========================================================================
%========================================================================
if nargin == 4
    flagtheoreticalSTDs = 0;
end  
%=== for theoretical STD computation
if flagtheoreticalSTDs
    allT.TUUonUR   = linspace(0.6,2,100);
    allT.TURonRR   = linspace(0.6,2,100);
    allT.MSC       = linspace(0.6,1,100);
    allT.phase     = linspace(-pi,pi,100);
    alphaCIforRsup = (1-normcdf(1))*2; % ci at 1 sigma
end
%=== for SCP analysis
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
        Lfft, overlapDFT, ...
        ratioDFT2SCP,overlapSCP, Fs_Hz, windowingshape);
    alltimes_sec{ifilter} = time_temp;
    %================================================================
    [estimRsup, estimRinf, spectralmatrix, allMSCs, Nsupthreshold] ...
        = estimSUT(allSDs, MSCthreshold);
    SUTs(ifilter).spectralmatrix = spectralmatrix;
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
            nansum(SUTs(ifilter).allMSCs.indexcst(...
            SUTs(ifilter).indexinsidefreqband(1):...
            SUTs(ifilter).indexinsidefreqband(2),:),2);
    end
    indexfrqinside = ...
        SUTs(ifilter).indexinsidefreqband(1):SUTs(ifilter).indexinsidefreqband(2);
    Lindexfrqinside = length(indexfrqinside);
    %====== STDs computation using the spectral matrices
    if flagtheoreticalSTDs
        for indfreq_ii = 1:Lindexfrqinside
            ifreq_ii = indexfrqinside(indfreq_ii);
            elementofR_ii = spectralmatrix(:,ifreq_ii);
            if any(isnan(elementofR_ii))
                SUTs(ifilter).theoSTDmodonRsup(indfreq_ii) = NaN;
                SUTs(ifilter).theoSTDphase_rd(indfreq_ii) = NaN;
            else
                RR_ii = [elementofR_ii(1) elementofR_ii(3) ; ...
                    elementofR_ii(3)' elementofR_ii(2)];
                [statUUonUR, statURonRR, statMSC, theoSTDphase_rd] = ...
                    theoreticalStats(allT, RR_ii,ratioDFT2SCP, ...
                    alphaCIforRsup);
                SUTs(ifilter).theoSTDmodonRsup(indfreq_ii) = ...
                    diff(statUUonUR.CI)/2;
                SUTs(ifilter).theoSTDphase_rd(indfreq_ii) = theoSTDphase_rd;
            end
        end
    end
end
%========================================================================
function [sigout, filterbank] = fbank(sigin,filtercharact,Fs_Hz)
%========================================================================
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
%========================================================================
function [allSDs, time_sec, frqsFFT_Hz] = ...
    estimSCP(xU,xR,Lfft,overlapFFT, ...
    ratioDFT2SCP, overlapSD, Fs_Hz, smoothwindow)
%========================================================================
% Perform the spectral components of the two signals xU et xR.
%========================================================================
% The code uses the Welch's approach: on a given time window of signals
% DFT is performed. Sucessive time windows can overlapped.
% Then the specral components are obtained by averaging some 
% successive DFT blocks.
%
% Inputs:
%    xU: signal observed on the SUT (T x 1)
%    xR: signal observed on the SREF (T x 1)
%    Lfft: length of FFTs
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
%========================================================================
%========================================================================
xU   = xU(:);
xR   = xR(:);
N    = length(xU);
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
%========================================================================
% the normalisation below is
% not useful if we only consider PSD ratios
Hwin = Hwin *sqrt(Lfft/(Hwin'*Hwin));
%========================================================================
%   <-------- NaverageFFTs = 5 ------->
%  |******|******|******|******|******|
%      |******|******|******|******|
% Given NaverageFFTs and the duration to performing
% the spectral components, we derive the length
% of the FFT. Then we derive the total nb of FFT blocks
% in the full file, taking into account the overlap.
% NblocsFFT = Ttotal/shiftFFT
%       typically shiftFFT = 1/2
% We compute all FFTs before to averaging.
% In another langage it is not necessary
% to proceed in such a way.
%
%========================================================================
sqrtLfft    = sqrt(Lfft);
shiftSignal = fix((1-overlapFFT)*Lfft);
NblocksFFT  = fix((N-(Lfft-shiftSignal))/shiftSignal);
allFFTsRR   = zeros(Lfft,NblocksFFT);
allFFTsUU   = zeros(Lfft,NblocksFFT);

for ibF  = 1:NblocksFFT
    ibT  = (ibF-1)*shiftSignal+(1:Lfft);
    xU_i = xU(ibT) .* Hwin;
    xU_i = xU_i-mean(xU_i);
    xR_i = xR(ibT) .* Hwin;
    xR_i = xR_i-mean(xR_i);
    allFFTsUU(:,ibF) = fft(xU_i,Lfft)/sqrtLfft;
    allFFTsRR(:,ibF) = fft(xR_i,Lfft)/sqrtLfft;
end
NaverageFFT = fix(ratioDFT2SCP/(1-overlapFFT))-1;
NaverageFFT = max([NaverageFFT,ratioDFT2SCP]);
shiftFFTs   = fix((1-overlapSD)*NaverageFFT);
allSDs      = struct;
time_sec    = struct;
NshiftFFTs_with_overlap = max([shiftFFTs,1]); 
NSD       = fix(NblocksFFT/NshiftFFTs_with_overlap);
for ibB=1:NSD,
    indB1 = (ibB-1)*NshiftFFTs_with_overlap+1;
    indB2 = indB1+NaverageFFT-1;
    indB  = fix(indB1):fix(indB2);
    indB  = indB(indB<= NblocksFFT);
    allSDs(ibB).RR  = mean(abs(allFFTsRR(:,indB)) .^ 2,2);
    allSDs(ibB).UU  = mean(abs(allFFTsUU(:,indB)) .^ 2,2);
    allSDs(ibB).UR  = mean(allFFTsRR(:,indB) .* ...
        conj(allFFTsUU(:,indB)),2);
    allSDs(ibB).MSC = (abs(allSDs(ibB).UR) .^2) ./...
        (allSDs(ibB).RR .* allSDs(ibB).UU);
    allSDs(ibB).Rsup  = allSDs(ibB).UU ./ allSDs(ibB).UR;
    allSDs(ibB).Rinf  = allSDs(ibB).UR ./ allSDs(ibB).RR;
    allSDs(ibB).det   = (abs(allSDs(ibB).RR) .* abs(allSDs(ibB).UU)) ...
        -(abs(allSDs(ibB).UR) .^2);
end
frqsFFT_Hz = (0:Lfft-1)*Fs_Hz/Lfft;
time_sec.FFT = ((0:NblocksFFT-1)+1/2)*shiftSignal/Fs_Hz;
time_sec.SD =  ((0:NSD-1)+1/2)*(shiftSignal/Fs_Hz)* ...
    NshiftFFTs_with_overlap;
time_sec.signals =  ((0:N-1)+1/2)/Fs_Hz;
%========================================================================
function [Rsup, Rinf, SCPsupeta, MSC, Nsupthreshold] ...
    = estimSUT(allSDs, MSCThreshold)
%========================================================================
% synopsis:
% estimate the stastitics of interest from the spectral
%     components. Then extract the values over the
%     threshold.
%
% [Rsup, Rinf, SpMatrix, MSC, Nsupthreshold] ...
%    = estimSUT(allSDs, MSCThreshold)
%========================================================================
%==== You have to change the 'SD' in 'SCP' to be homogeous
%==== You have to chamge the 'overeta' by 'cst'
%========================================================================
% Inputs:
%        allSDs, also called SCP in other programs:
%        SD (Spectral Density) = SCP (Spectral ComPonent)
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
%            Rsup.modcst: modulus of Rsup estimated on Rsup.tabmod above
%            threshold
%
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
%           MSC.tabweightMSC: table of weights for averaging the Rsup (see formula
%           2.14 to 2.17 in calibreport.pdf)
%       
%       Nsupthreshold: number of values above the threshold 
%
%
%========================================================================
%========================================================================

%========================================================================
%================ if weightingflag = 1
% a ponderation is applied using the estimated variance
% of the MSC estimate
%========================================================================
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
%===
% the 3 elements of spectral matrix at each frequency 
% are performed by an averaging on the full time slots 
% where the MSC is above the threshold. A trimmed mean 
% or a weighted mean could be applied.
% The spectral matrices are useful to compute the
% theoretical STDs.
%
SCPsupeta(1,:) = nanmean(tabRRwithMSCsupeta,2);
SCPsupeta(2,:) = nanmean(tabUUwithMSCsupeta,2);
SCPsupeta(3,:) = nanmean(tabURwithMSCsupeta,2);

%========================================================================
%================= on Rsup, Rinf
%========================================================================
tabHUUUR    = [allSDs.Rsup];
tabHURRR    = [allSDs.Rinf];
%========================================================================
%================= on Rsup
%===== real/imaginary part without constraint on MSC
%========================================================================
% Remark: in Matlab we can also use:
% sqrt(mean(real(z))^2+mean(imag(z))^2) = abs(mean(z))
% and 
% atan2(mean(imag(z)),mean(real(z))) = angle(mean(z))
%========================================================================
tabrealHUUUR  = real(tabHUUUR);
tabimagHUUUR  = imag(tabHUUUR);
tabmodHUUUR   = sqrt(tabrealHUUUR .^2 + tabimagHUUUR .^2);
tabphaseHUUUR = atan2(tabimagHUUUR,tabrealHUUUR);

%========================================================================
%================= on Rinf
%===== real/imaginary part without constraint on MSC
%========================================================================
tabrealHURRR  = real(tabHURRR);
tabimagHURRR  = imag(tabHURRR);
tabmodHURRR   = sqrt(tabrealHURRR .^2 + tabimagHURRR .^2);
tabphaseHURRR = atan2(tabimagHURRR, tabrealHURRR);

%========================================================================
%================= Rsup HAT by averaging
%========================================================================
hatrealHUUUR = nanmean(tabrealHUUUR,2);
hatimagHUUUR = nanmean(tabimagHUUUR,2);
%===== mod/phase part without constraint
hatmodHUUUR   = sqrt(hatrealHUUUR .^2 + hatimagHUUUR .^2);
hatphaseHUUUR = atan2(hatimagHUUUR, hatrealHUUUR);

%========================================================================
%================= Rinf HAT by averaging
%========================================================================
hatrealHURRR = nanmean(tabrealHURRR,2);
hatimagHURRR = nanmean(tabimagHURRR,2);
%===== mod/phase part without constraint
hatmodHURRR   = sqrt(hatrealHURRR .^2 + hatimagHURRR .^2);
hatphaseHURRR = atan2(hatimagHURRR, hatrealHURRR);

%========================================================================
%================ Rsup with constraint on the MSC
%========================================================================
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
%===== phase with constraint
tabphaseHUUURwithMSCsupeta = atan2(tabimagHUUURwithMSCsupeta,...
    tabrealHUUURwithMSCsupeta);

%========================================================================
%================ Rinf with constraint on the MSC
%========================================================================
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
%===== phase with constraint
tabphaseHURRRwithMSCsupeta = atan2(tabimagHURRRwithMSCsupeta,...
    tabrealHURRRwithMSCsupeta);


%========================================================================
%======= useful to perform the weights
RsuptheowithMSCsupeta = tabmodHUUURwithMSCsupeta .* ...
    tabmodHURRRwithMSCsupeta;

%========================================================================
%================ if weightingflag = 1
% a ponderation is applied using the estimated variance
% of the MSC estimate
%========================================================================
tabweightMSCsupeta = (tabMSCwithMSCsupeta .^2) ./  ...
    (1-tabMSCwithMSCsupeta) .* RsuptheowithMSCsupeta;
weightMSCinfeta = 1 ./  ...
    (1-tabMSCwithMSCsupeta) .* RsuptheowithMSCsupeta;
%========================================================================
%========== Rsup HAT with constraint and weighting coeffs
%========================================================================
if weightingflag
    hatrealHUUURwithMSCsupeta = ...
        nansum(tabrealHUUURwithMSCsupeta .* tabweightMSCsupeta,2) ...
        ./ nansum((tabweightMSCsupeta),2);
    
    hatimagHUUURwithMSCsupeta = ...
        nansum(tabimagHUUURwithMSCsupeta .* tabweightMSCsupeta,2) ...
        ./ nansum((tabweightMSCsupeta),2);
else
    hatrealHUUURwithMSCsupeta = ...
        nanmean(tabrealHUUURwithMSCsupeta,2);
    
    hatimagHUUURwithMSCsupeta = ...
        nanmean(tabimagHUUURwithMSCsupeta,2);
end

%===== perform module and phase part with constraint
hatmodHUUURwithMSCsupeta = sqrt(...
    hatrealHUUURwithMSCsupeta .^2+...
    hatimagHUUURwithMSCsupeta .^2);

hatphaseHUUURwithMSCsupeta = atan2(...
    hatimagHUUURwithMSCsupeta,...
    hatrealHUUURwithMSCsupeta);

%===== perform STD on module and phase
stdmodHUUURwithMSCsupeta   = nanstd(tabmodHUUURwithMSCsupeta,[],2);
stdphaseHUUURwithMSCsupeta = nanstd(tabphaseHUUURwithMSCsupeta,[],2);

%========================================================================
%========== Rinf HAT with constraint and weighting coeffs
%========================================================================
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

%===== perform module and phase part with constraint
hatmodHURRRwithMSCsupeta = sqrt(...
    hatrealHURRRwithMSCsupeta .^2+...
    hatimagHURRRwithMSCsupeta .^2);

hatphaseHURRRwithMSCsupeta = atan2(...
    hatimagHURRRwithMSCsupeta,...
    hatrealHURRRwithMSCsupeta);

%===== perform STD on module and phase
stdmodHURRRwithMSCsupeta   = nanstd(tabmodHURRRwithMSCsupeta,[],2);
stdphaseHURRRwithMSCsupeta = nanstd(tabphaseHURRRwithMSCsupeta,[],2);

%========================================================================
%=========================== SAVE =======================================
Rsup.tabmod      = tabmodHUUUR;
Rsup.tabmodcst   = tabmodHUUURwithMSCsupeta;
Rsup.tabphase    = tabphaseHUUUR;
Rsup.tabphasecst = tabphaseHUUURwithMSCsupeta;

Rsup.mod         = hatmodHUUUR;
Rsup.phase       = hatphaseHUUUR;
Rsup.modcst      = hatmodHUUURwithMSCsupeta;
Rsup.stdmodcst   = stdmodHUUURwithMSCsupeta;
Rsup.phasecst    = hatphaseHUUURwithMSCsupeta;
Rsup.stdphasecst = stdphaseHUUURwithMSCsupeta;

Rinf.tabmod      = tabmodHURRR;
Rinf.tabmodcst   = tabmodHURRRwithMSCsupeta;
Rinf.tabphase    = tabphaseHURRR;
Rinf.tabphasecst = tabphaseHURRRwithMSCsupeta;

Rinf.mod         = hatmodHURRR;
Rinf.phase       = hatphaseHURRR;
Rinf.modcst      = hatmodHURRRwithMSCsupeta;
Rinf.stdmodcst   = stdmodHURRRwithMSCsupeta;
Rinf.phasecst    = hatphaseHURRRwithMSCsupeta;
Rinf.stdphasecst = stdphaseHURRRwithMSCsupeta;

MSC.tab          = tabMSC;
MSC.tabcst       = tabMSCwithMSCsupeta;
MSC.indexcst     = indextabMSCsupthreshold;
MSC.tabweightMSC    = tabweightMSCsupeta;
%========================= END of analysis ===============================
%=========================================================================
function [STATSmodUUonRU, STATSmodURonRR, statMSC, ...
    STDphase_rd, statPhasepsd] = ...
    theoreticalStats(allT, Smatrix, N, alphaCI)
%=========================================================================
% Compute, for a given frequency, the probability density functions, 
% the STDs and the CI (at alphaCI) of:
%    - the ratio modules:
%
%               |Smatrix(1,2)|                  Smatrix(1,1)
%      Rinf =  ---------------   and   Rsup =  --------------
%                Smatrix(2,2)                  |Smatrix(2,1)|
%    - the phase
%                phase = arg Smatrix(1,2)
%    - the MSC 
%                          |Smatrix(1,2)|^2
%                MSC = -------------------------
%                       Smatrix(1,1)Smatrix(2,2)
% Rk:  |Smatrix(1,2)| = |Smatrix(2,1)| because Smatrix is positive
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
    exp(-(allT.phase-phase_HUminusHR_rad) .^2/...
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





