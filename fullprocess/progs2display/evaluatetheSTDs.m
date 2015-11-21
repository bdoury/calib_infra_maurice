% this program performs theoretical STDs for a group of randomly
% selected pairs of days, and also the STDs obtained on the 
% estimates.
%
%=========================================================================
clear
addpath ../ZZtoolbox/
directorysignals    = '../../../../AAdataI26calib/';

%==== this directory contains the parameters evalauted by the
% program estimationwithFB.m
directoryinputresults = '../AAresultswithFB98bis/';

%======== to compute theoretical integrals
allT.TUUonUR   = linspace(0.6,2,100);
allT.TURonRR   = linspace(0.6,2,100);
allT.MSC       = linspace(0.6,1,100);
allT.phase     = linspace(-pi,pi,100);
alphaSTDforRsup = (1-normcdf(1))*2; % ci at 1 sigma

randomlydoubledaynumber = 15;

for ihc = 2
    numfig = ihc;
    % list of the files from 1 to nbmats
    % if you want a name type fileswithdotmat(#)
    fileswithdotmat = dir(sprintf('../%ss%i/s%iy*.mat',...
        directorysignals,ihc,ihc));
    comload = sprintf('load %sresultssta26sensor%i',directoryinputresults,ihc);
    eval(comload);
    switch ihc
        case 1
            remainindex = [1:34 36:61  63:70]; % 2015/10/13
        case 4
            remainindex = [1:nbmats]; % 2015/10/11
        otherwise
            remainindex = [1:nbmats];
    end
    %%
    
    allRatioSupPfilters          = allRatioSupPfilters(:,remainindex);
    allSTDmodRatioSupPfilters    = allSTDmodRatioSupPfilters(:,remainindex);
    allSTDphaseRatioSupPfilters  = allSTDphaseRatioSupPfilters(:,remainindex);
    
    allRatioInfPfilters          = allRatioInfPfilters(:,remainindex);
    allSTDmodRatioInfPfilters    = allSTDmodRatioInfPfilters(:,remainindex);
    allSTDphaseRatioInfPfilters  = allSTDphaseRatioInfPfilters(:,remainindex);
    allmeanMSCcstPfilters        = allmeanMSCcstPfilters(:,remainindex);
    nbofvaluesoverthreshold      = nbofvaluesoverthreshold(:,remainindex);
    allScpPfilters               = allScpPfilters(:,:,remainindex);
    %%
    permutenbmats = randperm(length(remainindex));
    indrandomlychosen = permutenbmats(1:randomlydoubledaynumber);
    allRatioSupPfilters = ...
        allRatioSupPfilters(:,indrandomlychosen);
    allSTDmodRatioSupPfilters = ...
        allSTDmodRatioSupPfilters(:,indrandomlychosen);
    allSTDphaseRatioSupPfilters = ...
        allSTDphaseRatioSupPfilters(:,indrandomlychosen);
    allmeanMSCcstPfilters = ...
        allmeanMSCcstPfilters(:,indrandomlychosen);
    nbofvaluesoverthreshold = ...
        nbofvaluesoverthreshold(:,indrandomlychosen);
    allScpPfilters = allScpPfilters(:,:,indrandomlychosen);
    
    STDmodRatioPfilters_ave        = nanmedian(allSTDmodRatioSupPfilters,2);
    STDphaseRatioPfilters_ave      = nanmedian(allSTDphaseRatioSupPfilters,2);
    
    %== sort in increasing order
    [allfrqsPfiltersS, inds]       = sort(allfrqsPfilters);
    allRatioPfiltersS              = allRatioSupPfilters(inds,:);
    allmeanMSCcstPfiltersS         = allmeanMSCcstPfilters(inds,:);
    STDmodRatioPfilters_aveS       = STDmodRatioPfilters_ave(inds);
    STDphaseRatioPfilters_aveS     = STDphaseRatioPfilters_ave(inds);
    nbofvaluesoverthresholdS       = nbofvaluesoverthreshold(inds,:);
    allScpPfiltersS                = allScpPfilters(:,inds,:);
    
    %== unique frequency value
    [allfrqsPfiltersUS, inda]      = unique(allfrqsPfiltersS);
    allRatioPfiltersUS             = allRatioPfiltersS(inda,:);
    allmeanMSCcstPfiltersUS        = allmeanMSCcstPfiltersS(inda,:);
    STDmodRatioPfilters_aveUS      = STDmodRatioPfilters_aveS(inda);
    STDphaseRatioPfilters_aveUS    = STDphaseRatioPfilters_aveS(inda);
    nbofvaluesoverthresholdUS      = nbofvaluesoverthresholdS(inda,:);
    allScpPfiltersUS               = allScpPfiltersS(:,inda,:);
    
    %== without 0 frequency values
    indz = find(not(allfrqsPfiltersUS==0));
    allfrqsPfiltersUSZ             = allfrqsPfiltersUS(indz);
    RatioPfiltersUSZ               = allRatioPfiltersUS(indz,:);
    allmeanMSCcstPfiltersUSZ       = allmeanMSCcstPfiltersUS(indz,:);
    STDmodRatioPfilters_aveUSZ     = STDmodRatioPfilters_aveUS(indz);
    STDphaseRatioPfilters_aveUSZ   = STDphaseRatioPfilters_aveUS(indz);
    nbofvaluesoverthresholdUSZ     = nbofvaluesoverthresholdUS(indz,:);
    allScpPfiltersUSZ              = allScpPfiltersS(:,indz,:);
    
    %====== absolute and arg of the ratios
    modRatioPfiltersUSZ            = abs(RatioPfiltersUSZ);
    phaseRatioPfiltersUSZ_rd       = angle(RatioPfiltersUSZ);
    %====== averaging by MEDIAN to avoid outliers
    trimmeanmodRatioPfiltersUSZ    = trimmean(modRatioPfiltersUSZ,30,2);
    meanmodRatioPfiltersUSZ        = nanmean(modRatioPfiltersUSZ,2);
    meanphasePfiltersUSZ_rd        = trimmean(phaseRatioPfiltersUSZ_rd,30,2);
    %====== STDs
    STDmodPfiltersUSZ              = nanstd(modRatioPfiltersUSZ,[],2);
    STDphasePfiltersUSZ_rd         = nanstd(phaseRatioPfiltersUSZ_rd,[],2);
    
    ICallRatioPfiltersUSZ          = STDmodPfiltersUSZ ./ ...
        sqrt(sum(nbofvaluesoverthresholdUSZ,2));
    
    ICallRatioPfiltersUSZbis       = STDmodRatioPfilters_aveUSZ ./ ...
        sqrt(sum(nbofvaluesoverthresholdUSZ,2));
end
%%
NaverageFFTs   = filtercharact(1).ratioDFT2SCP;
LfqsUSZ       = length(allfrqsPfiltersUSZ);
twolistsSTDs  = zeros(LfqsUSZ,2);
STDmodtheoUSZ    = zeros(LfqsUSZ,1);
STDphasetheoUSZ  = zeros(LfqsUSZ,1);

for indfq=1:LfqsUSZ
    SCP_ifq = nanmean(allScpPfiltersUSZ(:,indfq,:),3);
    if any(isnan(SCP_ifq))
        STDmodtheoUSZ(indfq)=NaN;
    else
        RR_ifq = [SCP_ifq(1) SCP_ifq(3)';SCP_ifq(3) SCP_ifq(2)];
        [statUUonUR, statURonRR, statMSC, stdPhase_degree]    = ...
            statsRatiosHbis(allT, RR_ifq, NaverageFFTs, alphaSTDforRsup);
        STDmodtheoUSZ(indfq) = diff(statUUonUR.CI)/2;
        STDphasetheoUSZ(indfq) = stdPhase_degree;
    end
end

% [nanmean(STDmodtheo), nanstd(STDmodtheo)],
% [nanmean(STDmodRatioPfilters_aveUSZ), nanstd(STDmodRatioPfilters_aveUSZ)],
% [nanmean(STDmodPfiltersUSZ), nanstd(STDmodPfiltersUSZ)]
%%
% probaIC 70%, 1.6449 for 90%  and 1.96 for 95%
probaIC = 0.95;
alphaCI = -norminv((1-probaIC)/2);

sumnbofvaluesoverthresholdUSZ = sum(nbofvaluesoverthresholdUSZ,2);
CItheoUSZ = alphaCI * STDmodtheoUSZ ./ sumnbofvaluesoverthresholdUSZ;
CIPfiltersUSZ = alphaCI * STDmodPfiltersUSZ ./ sumnbofvaluesoverthresholdUSZ;
CIRatioPfilters_aveUSZ = alphaCI * STDmodRatioPfilters_aveUSZ ./ ...
    sumnbofvaluesoverthresholdUSZ;
CISTDphasetheoUSZ = STDphasetheoUSZ ./ sumnbofvaluesoverthresholdUSZ;
figure(2)
subplot(311)
loglog(allfrqsPfiltersUSZ,sumnbofvaluesoverthresholdUSZ,'.')
grid on
subplot(312)
loglog(allfrqsPfiltersUSZ, [CItheoUSZ ...
    CIPfiltersUSZ CIRatioPfilters_aveUSZ],'.')
grid on
subplot(313)
loglog(allfrqsPfiltersUSZ, CISTDphasetheoUSZ,'.')
grid on
%==================================
%
% plot(allT.TUUonUR,statURonRR.pdf,'.-')
% hold on
% plot(allT.TUUonUR,statUUonUR.pdf,'.-r')
% hold off

