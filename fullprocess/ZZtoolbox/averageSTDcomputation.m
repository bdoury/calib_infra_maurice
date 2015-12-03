function [STDmodtheo, STDphasetheo_rd] = ...
    averageSTDcomputation(allT, M, SPs, ...
    allfrqsPfilters, alphaIC)
%=========================================================================
% performs the theoretical STDs of the Sup ratio module and
% the Sup ratio phase.
% Inputs:
%    - allT: range of values to compute theoretical integrals,
%      typical ranges:
%      allT.TUUonUR    = linspace(0.6,2,100);
%      allT.TURonRR    = linspace(0.6,2,100);
%      allT.MSC        = linspace(0.6,1,100);
%      allT.phase_rd   = linspace(-pi,pi,100);
%    - M: window length for averaging FFT frames
%    - SPs : spectral matrix components
%            3 x K, where 
%           SPs(1,x) = autospectrum on 1
%           SPs(2,x) = autospectrum on 2
%           SPs(3,x) = crossspectrum between 1 and 2
%    - alphaIC: level of confidence, typically 5%
%
% Used functions: theoreticalStats.m
%                 norminv.m from Matlab
%=========================================================================
% M      = filtercharact(1).ratioDFT2SCP;
% meanSPs = squeeze(nanmean(allScpPfilters,3));
Lfqs             = length(allfrqsPfilters);
STDmodtheo       = NaN(Lfqs,1);
STDphasetheo_rd  = NaN(Lfqs,1);

for indfq=1:Lfqs
    SCP_ifq = SPs(:,indfq);
    if any(isnan(SCP_ifq))
        STDmodtheo(indfq)=NaN;
    else
        RR_ifq = [SCP_ifq(1) SCP_ifq(3)';SCP_ifq(3) SCP_ifq(2)];
        [statUUonUR, statURonRR, statMSC, stdPhase_rd]    = ...
            theoreticalStats(allT, RR_ifq, M, alphaIC);
        STDmodtheo(indfq) = diff(statUUonUR.CI)/2;
        STDphasetheo_rd(indfq) = stdPhase_rd;
    end
end