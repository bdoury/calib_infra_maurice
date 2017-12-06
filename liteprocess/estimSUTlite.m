function [Rsup, freqslin, STDmoduleRlin, ...
    STDphaseRlin_rd, nboverTHlin] = ...
    estimSUTlite ...
    (signals, structfiltercharacteristics, frequencylist_Hz, ...
    Fs_Hz, MSCthreshold, trimpercent)
%========================================================================
% Synopsis:
% [Rsup, freqslin, STDmoduleRlin, ...
%     STDphaseRlin_rd, nboverTHlin] = ...
%     estimSUTlite ...
%     (signals, structfiltercharacteristics, frequencylist_Hz, ...
%     Fs_Hz, MSCthreshold, trimpercent)
%===============
% Inputs:
%     - signals : T x 2
%     - structfiltercharacteristics (FB structure)
%           see document
%     - frequencylist_Hz: array N x 1 of the selected frequencies
%       in Hz. N can take any value under Fs_Hz/2
%     - Fs_Hz: sampling frequency in Hz
%     - MSCthreshold:
%     - trimpercent: percent of values keptfor averaging
%===============
% Outputs:
%     - Rsup: array N x 1 of the estimated ratios
%     - freqslin: array N x 1 of the selected frequencies
%       in Hz. almost the same as frequencylist_Hz, except if some 
%       are outside of the FB bandwidths.
%     - STDmoduleR: array N x 1 of the STD on the module of Rsup
%     - STDphaseR_rd: array N x 1 of the STD on the phase of Rsup
%     - nboverTH: array N x 1 of the number of values over the threshold
%========================================================================

nbfrequencies          = length(frequencylist_Hz);
Pfilter                = length(structfiltercharacteristics);
frequenciesinfilter_Hz = cell(Pfilter,1);
nbfreqsbyfilter        = NaN(Pfilter,1);

%=== determine the frequencies inside the bank filters
%  in such a way that all frequencies are only in
%  ONE filter band
frequencylist_Hz_ii = frequencylist_Hz;
nbfrequencies_ii    = nbfrequencies;
for idfilter=1:Pfilter
    fqlow_Hz    = structfiltercharacteristics(idfilter).Wlow_Hz;
    fqhigh_Hz   = structfiltercharacteristics(idfilter).Whigh_Hz;
    cp=0;
    for idf=1:nbfrequencies_ii
        if and(frequencylist_Hz_ii(idf)>fqlow_Hz, ...
                frequencylist_Hz_ii(idf)<=fqhigh_Hz)
            cp=cp+1;
            frequenciesinfilter_Hz{idfilter}(cp) = ...
                frequencylist_Hz_ii(idf);
        end
    end
    nbfreqsbyfilter(idfilter) = cp;
    frequencylist_Hz_ii = ...
        setdiff(frequencylist_Hz_ii,frequenciesinfilter_Hz{idfilter});
    nbfrequencies_ii = length(frequencylist_Hz_ii);
end
nbofallfrequencies = sum(nbfreqsbyfilter);
%========== we perform the filter coeffiecient from the structure
%           denoted structfiltercharacteristics
% using the Matlab functions as BUTTER.M
filterbankcoeff = cell(Pfilter,1);
for ifilter = 1:Pfilter
    fname   = structfiltercharacteristics(ifilter).designname;
    forder  = structfiltercharacteristics(ifilter).Norder;
    fqlow   = structfiltercharacteristics(idfilter).Wlow_Hz/Fs_Hz;
    fqhigh  = structfiltercharacteristics(idfilter).Whigh_Hz/Fs_Hz;
    switch fname
        case 'fir1'
            fdesign = sprintf('filnum = %s(%i,[%5.8f,%5.8f]);',...
                fname,forder,2*fqlow,2*fqhigh);
            filden = 1;
        case 'butter'
            fdesign = sprintf('[filnum,filden] = %s(%i,[%5.8f %5.8f]);',...
                fname,forder,2*fqlow,2*fqhigh);
        case 'cheby1'
            fdesign = sprintf('[filnum,filden] = %s(%i,%i,[%5.8f %5.8f]);',...
                fname,forder,0.02,2*fqlow,2*fqhigh);
    end
    eval(fdesign)
    filterbankcoeff{ifilter}.num = filnum;
    filterbankcoeff{ifilter}.den = filden;
end
%========== we perform the shape window from the structure
%           denoted structfiltercharacteristics
% using the Matlab functions as HANN.M
windshape = cell(Pfilter,1);
for ifilter = 1:Pfilter
    windowshapename = structfiltercharacteristics(ifilter).windowshape;
    SCPperiod_sec   = structfiltercharacteristics(ifilter).SCPperiod_sec;
    ratioDFT2SCP    = structfiltercharacteristics(ifilter).ratioDFT2SCP;
    lengthDFT       = fix(SCPperiod_sec*Fs_Hz/ratioDFT2SCP);
    switch windowshapename
        case 'hann'
            windshape{ifilter} = hann(lengthDFT,'periodic');
            windshape{ifilter} = windshape{ifilter} / ...
                sqrt(sum(windshape{ifilter} .^2));
    end
end
%==== pre-computation of the exponentials used by
% the direct DFTs
EXPV                = cell(Pfilter,1);
for ifilter = 1:Pfilter
    SCPperiod_sec   = structfiltercharacteristics(ifilter).SCPperiod_sec;
    ratioDFT2SCP    = structfiltercharacteristics(ifilter).ratioDFT2SCP;
    lengthDFT       = fix(SCPperiod_sec*Fs_Hz/ratioDFT2SCP);
    DFTindex        = (0:lengthDFT-1)'/Fs_Hz;
    EXPV{ifilter}   = exp(-2j*pi*DFTindex*frequenciesinfilter_Hz{ifilter});
end
%============================================
%============================================
Nsignals   = size(signals,1);
R          = cell(Pfilter,1);
STDmoduleR = cell(Pfilter,1);
STDphaseR  = cell(Pfilter,1);
nboverTH   = cell(Pfilter,1);
for ifilter = 1:Pfilter
    filnum = filterbankcoeff{ifilter}.num;
    filden = filterbankcoeff{ifilter}.den;
    filteredsignals = filter(filnum,filden,signals);
    SCPperiod_sec   = structfiltercharacteristics(ifilter).SCPperiod_sec;
    ratioDFT2SCP    = structfiltercharacteristics(ifilter).ratioDFT2SCP;
    overlapDFT      = structfiltercharacteristics(ifilter).overlapDFT;
    % Computation
    lengthDFT       = fix(SCPperiod_sec*Fs_Hz/ratioDFT2SCP);
    lengthSCP       = fix(SCPperiod_sec*Fs_Hz);
    DFTshift        = fix((1-overlapDFT)*lengthDFT);
    NSCPwindows     = fix(Nsignals/Fs_Hz/SCPperiod_sec);
    sigauxW         = zeros(lengthDFT,2);
    
    SCP_ifreq11     = zeros(nbfreqsbyfilter(ifilter),NSCPwindows-1);
    SCP_ifreq22     = zeros(nbfreqsbyfilter(ifilter),NSCPwindows-1);
    SCP_ifreq12     = zeros(nbfreqsbyfilter(ifilter),NSCPwindows-1);
    
    for iwindowSCP  = 1:NSCPwindows-1
        id0   = (iwindowSCP-1)*lengthSCP;
        id1   = 0;
        cpDFT = 0;
        while id1<id0+lengthSCP-lengthDFT
            cpDFT  = cpDFT+1;
            id1    = id0 + (cpDFT-1)*DFTshift+1;
            id2    = id1+lengthDFT-1;
            sigaux = filteredsignals(id1:id2,:);
            sigauxW(:,1) = sigaux(:,1) .* windshape{ifilter};
            sigauxW(:,2) = sigaux(:,2) .* windshape{ifilter};
            for ifreq = 1:nbfreqsbyfilter(ifilter)
                X_ifreq1 = sum(sigauxW(:,1) .* EXPV{ifilter}(:,ifreq));
                X_ifreq2 = sum(sigauxW(:,2) .* EXPV{ifilter}(:,ifreq));
                SCP_ifreq11(ifreq,iwindowSCP) = SCP_ifreq11(ifreq,iwindowSCP) + ...
                    X_ifreq1 .* conj(X_ifreq1);
                SCP_ifreq22(ifreq,iwindowSCP) = SCP_ifreq22(ifreq,iwindowSCP) + ...
                    X_ifreq2 .* conj(X_ifreq2);
                SCP_ifreq12(ifreq,iwindowSCP) = SCP_ifreq12(ifreq,iwindowSCP) + ...
                    (X_ifreq1) .* conj(X_ifreq2);
            end
        end
    end
    tabMSC_ifilter   = (abs(SCP_ifreq12) .^2) ./ ...
        (SCP_ifreq11 .* SCP_ifreq22);
    ind_ifilter_cst = (tabMSC_ifilter>MSCthreshold);
    tabMSC_ifilter_cst = NaN(size(tabMSC_ifilter));
    tabMSC_ifilter_cst(ind_ifilter_cst) = ...
        tabMSC_ifilter(ind_ifilter_cst);
    
    tabRsup_ifilter  = SCP_ifreq11 ./ conj(SCP_ifreq12);
    tabRsup_ifilter_cst = ...
        NaN(size(tabRsup_ifilter))+1j*NaN(size(tabRsup_ifilter));
    tabRsup_ifilter_cst(ind_ifilter_cst) = ...
        tabRsup_ifilter(ind_ifilter_cst);
    tabRsup_ifilter_cst_trim = ...
        trimmeancomplex(tabRsup_ifilter_cst,trimpercent);    
    
    SCP_ifreq11_cst   = NaN(size(SCP_ifreq11));
    SCP_ifreq11_cst(ind_ifilter_cst) = SCP_ifreq11(ind_ifilter_cst);
    SCP_ifreq22_cst   = NaN(size(SCP_ifreq22));
    SCP_ifreq22_cst(ind_ifilter_cst) = SCP_ifreq22(ind_ifilter_cst);
    tabR1122_cst      = SCP_ifreq11_cst ./ SCP_ifreq22_cst;

    weightMSCsupeta  = (tabMSC_ifilter_cst .^2) ./  ...
        (1-tabMSC_ifilter_cst) .* tabR1122_cst;
    
    R_filter         = ...
        nansum(tabRsup_ifilter_cst_trim .* weightMSCsupeta,2) ...
        ./ nansum(weightMSCsupeta,2);

    R{ifilter}       = R_filter;
    nboverTH_ii      = nansum(ind_ifilter_cst,2);
    %===== perform STD on module and phase
    STDmoduleR{ifilter}   = nanstd(abs(tabRsup_ifilter_cst),[],2) ./...
        sqrt(nboverTH_ii);
    STDphaseR{ifilter}    = nanstd(angle(tabRsup_ifilter_cst),[],2) ./...
        sqrt(nboverTH_ii);
    nboverTH{ifilter}     = nboverTH_ii;
end
freqslin      = zeros(nbofallfrequencies,1);
Rsup          = zeros(nbofallfrequencies,1);
id2           = 0;
STDmoduleRlin = zeros(nbofallfrequencies,1);
STDphaseRlin_rd  = zeros(nbofallfrequencies,1);
nboverTHlin   = zeros(nbofallfrequencies,1);
for ifilter=1:Pfilter
    id1 = id2+1;
    id2 = id1+nbfreqsbyfilter(ifilter)-1;
    Rsup(id1:id2)     = R{ifilter};
    freqslin(id1:id2) = frequenciesinfilter_Hz{ifilter}';
    STDmoduleRlin(id1:id2) = STDmoduleR{ifilter};
    STDphaseRlin_rd(id1:id2) = STDphaseR{ifilter};
    nboverTHlin(id1:id2) =  nboverTH{ifilter};
end

%===================================================
%===================================================
%===================================================
function trimmedz = trimmeancomplex(z,trimpercent)

[ra,co] = size(z);
trimmedz = nan(ra,co);

for ira=1:ra 
    indout   = quadform(z(ira,:),trimpercent);
    trimmedz(ira,indout==1) = z(ira,indout==1);
end
%===================================================
function indout = quadform(z, apercent)
%===================================================
c     = -2*log(1-apercent);
z     = z(:);
N     = length(z);
meanz = nanmean(z);
zc    = z-ones(N,1)*meanz;
HH    = [real(z) imag(z)];
R     = nancov(HH);
rizc  = [real(zc), imag(zc)];
if or(sum(any(isnan(R)))>0,sum(any(isinf(R)))>0)
    indout = zeros(N,1);
elseif rank(R)==2
    Fm1   = sqrtm(R);
    % valp   = eig(R);
    % area   = sqrt(prod(valp))*c*pi;
    indout = zeros(N,1);
    for ii=1:N
        indout(ii) = rizc(ii,:) * Fm1 *rizc(ii,:)'<c;
    end
else
    indout = zeros(N,1);
end
%=====================================================


