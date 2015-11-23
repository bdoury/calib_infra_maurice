%====================== evaluatetheSTDs.m ================================
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
allT.TUUonUR    = linspace(0.6,2,100);
allT.TURonRR    = linspace(0.6,2,100);
allT.MSC        = linspace(0.6,1,100);
allT.phase_rd   = linspace(-pi,pi,100);
% critical region at +/-1 sigma
alphaSTDforRsup = (1-normcdf(1))*2; 
probaIC         = 0.90;
alphaCI         = -norminv((1-probaIC)/2);
drawnumber      = 15;

for ihc = 1
    numfig = ihc;
    % list of the files from 1 to nbmats
    % if you want a name type fileswithdotmat(#)
    fileswithdotmat = dir(sprintf('../%ss%i/s%iy*.mat',...
        directorysignals,ihc,ihc));
    comload = sprintf('load %sresultssta26sensor%i',directoryinputresults,ihc);
    eval(comload);
    switch ihc
        case 1
            remainindex = [1:34 36:61  63:70]; % 2015/08/09 and 2015/10/13
        case 4
            remainindex = [1:nbmats]; % 2015/10/11
        otherwise
            remainindex = [1:nbmats];
    end
    %%
    fileswithdotmatremain        = fileswithdotmat(remainindex);
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
    permutenbmats                = randperm(length(remainindex));
    indrandomlychosen            = permutenbmats(1:drawnumber);
    allRatioSupPfilters          = ...
        allRatioSupPfilters(:,indrandomlychosen);
    allSTDmodRatioSupPfilters    = ...
        allSTDmodRatioSupPfilters(:,indrandomlychosen);
    allSTDphaseRatioSupPfilters  = ...
        allSTDphaseRatioSupPfilters(:,indrandomlychosen);
    allmeanMSCcstPfilters        = ...
        allmeanMSCcstPfilters(:,indrandomlychosen);
    nbofvaluesoverthreshold      = ...
        nbofvaluesoverthreshold(:,indrandomlychosen);
    allScpPfilters               = ...
        allScpPfilters(:,:,indrandomlychosen);
    
    STDmodRatioPfilters_ave      = nanmedian(allSTDmodRatioSupPfilters,2);
    STDphaseRatioPfilters_ave    = nanmedian(allSTDphaseRatioSupPfilters,2);
    
    %== sort in increasing order
    [allfrqsPfiltersS, inds]     = sort(allfrqsPfilters);
    allRatioPfiltersS            = allRatioSupPfilters(inds,:);
    allmeanMSCcstPfiltersS       = allmeanMSCcstPfilters(inds,:);
    STDmodRatioPfilters_aveS     = STDmodRatioPfilters_ave(inds);
    STDphaseRatioPfilters_aveS   = STDphaseRatioPfilters_ave(inds);
    nbofvaluesoverthresholdS     = nbofvaluesoverthreshold(inds,:);
    allScpPfiltersS              = allScpPfilters(:,inds,:);
    
    %== unique frequency value
    [allfrqsPfiltersUS, inda]    = unique(allfrqsPfiltersS);
    allRatioPfiltersUS           = allRatioPfiltersS(inda,:);
    allmeanMSCcstPfiltersUS      = allmeanMSCcstPfiltersS(inda,:);
    STDmodRatioPfilters_aveUS    = STDmodRatioPfilters_aveS(inda);
    STDphaseRatioPfilters_aveUS  = STDphaseRatioPfilters_aveS(inda);
    nbofvaluesoverthresholdUS    = nbofvaluesoverthresholdS(inda,:);
    allScpPfiltersUS             = allScpPfiltersS(:,inda,:);
  
    %== without 0 frequency values
    indz                         = find(not(allfrqsPfiltersUS==0));
    allfrqsPfiltersUSZ           = allfrqsPfiltersUS(indz);
    RatioPfiltersUSZ             = allRatioPfiltersUS(indz,:);
    allmeanMSCcstPfiltersUSZ     = allmeanMSCcstPfiltersUS(indz,:);
    STDmodRatioPfilters_aveUSZ   = STDmodRatioPfilters_aveUS(indz);
    STDphaseRatioPfilters_aveUSZ = STDphaseRatioPfilters_aveUS(indz);
    nbofvaluesoverthresholdUSZ   = nbofvaluesoverthresholdUS(indz,:);
    allScpPfiltersUSZ            = allScpPfiltersS(:,indz,:);
    
    %====== absolute and arg of the ratios
    modRatioPfiltersUSZ          = abs(RatioPfiltersUSZ);
    phaseRatioPfiltersUSZ_rd     = angle(RatioPfiltersUSZ);
    %====== averaging by TRIMMEAN to avoid outliers
    trimmeanmodRatioPfiltersUSZ  = trimmean(modRatioPfiltersUSZ,30,2);
    meanmodRatioPfiltersUSZ      = nanmean(modRatioPfiltersUSZ,2);
    meanphasePfiltersUSZ_rd      = trimmean(phaseRatioPfiltersUSZ_rd,30,2);
    %====== STDs
    STDmodPfiltersUSZ            = nanstd(modRatioPfiltersUSZ,[],2);
    STDphasePfiltersUSZ_rd       = nanstd(phaseRatioPfiltersUSZ_rd,[],2);
    
    ICallRatioPfiltersUSZ        = STDmodPfiltersUSZ ./ ...
        sqrt(sum(nbofvaluesoverthresholdUSZ,2));
    
    ICallRatioPfiltersUSZbis     = STDmodRatioPfilters_aveUSZ ./ ...
        sqrt(sum(nbofvaluesoverthresholdUSZ,2));
end
%%
% this part can be removed because these values are
% performed on the last version of fbanalysis.m 
%====
NaverageFFTs     = filtercharact(1).ratioDFT2SCP;
LfqsUSZ          = length(allfrqsPfiltersUSZ);
twolistsSTDs     = zeros(LfqsUSZ,2);
STDmodtheoUSZ_rd    = zeros(LfqsUSZ,1);
STDphasetheoUSZ_rd  = zeros(LfqsUSZ,1);

for indfq=1:LfqsUSZ
    SCP_ifq = nanmean(allScpPfiltersUSZ(:,indfq,:),3);
    if any(isnan(SCP_ifq))
        STDmodtheoUSZ_rd(indfq)=NaN;
    else
        RR_ifq = [SCP_ifq(1) SCP_ifq(3)';SCP_ifq(3) SCP_ifq(2)];
        [statUUonUR, statURonRR, statMSC, stdPhase_rd]    = ...
            theoreticalStats(allT, RR_ifq, NaverageFFTs, alphaSTDforRsup);
        STDmodtheoUSZ_rd(indfq) = diff(statUUonUR.CI)/2;
        STDphasetheoUSZ_rd(indfq) = stdPhase_rd;
    end
end
sumnbofvaluesoverthresholdUSZ = sum(nbofvaluesoverthresholdUSZ,2);
CItheoUSZ = alphaCI * STDmodtheoUSZ_rd ./ sumnbofvaluesoverthresholdUSZ;
CIPfiltersUSZ = alphaCI * STDmodPfiltersUSZ ./ sumnbofvaluesoverthresholdUSZ;
CImodPfilters_aveUSZ = alphaCI * STDmodRatioPfilters_aveUSZ ./ ...
    sumnbofvaluesoverthresholdUSZ;

CISTDphasePfiltersUSZ_degree       = (180/pi) * alphaCI * STDphasePfiltersUSZ_rd./ ...
    sumnbofvaluesoverthresholdUSZ;
CISTDphasePfilters_aveUSZ_degree   = (180/pi) * alphaCI * STDphaseRatioPfilters_aveUSZ ./ ...
    sumnbofvaluesoverthresholdUSZ;
CISTDphasetheoUSZ_degree = (180/pi) * alphaCI * STDphasetheoUSZ_rd ./ sumnbofvaluesoverthresholdUSZ;
%%
figure(2)
subplot(311)
loglog(allfrqsPfiltersUSZ,sumnbofvaluesoverthresholdUSZ,'ob','markerfacec','b')
title(sprintf('IS26 -  sensor H%i, MSC threshold = %4.2f\nday number = %i',...
    ihc, MSCthreshold, 2*drawnumber),'fontname','times','fontsize',14)
ylabel('counts above the threshold','fontname','times','fontsize',12)
grid on
subplot(312)
loglog(allfrqsPfiltersUSZ, CItheoUSZ, 'ob','markerfacec','b')
hold on
loglog(allfrqsPfiltersUSZ, CImodPfilters_aveUSZ,'or','markerfacec','r')
hold off
grid on
ylabel(sprintf('CI at %i%s\non the module',probaIC*100,'%'), ...
    'fontname','times','fontsize',12)
hl=legend('theoretical','emprirically');
set(hl,'fontname','times','fontsize',12)

subplot(313)
loglog(allfrqsPfiltersUSZ, CISTDphasetheoUSZ_degree,'ob','markerfacec','b')
hold on
loglog(allfrqsPfiltersUSZ, CISTDphasePfilters_aveUSZ_degree,'or','markerfacec','r')
hold on
loglog(allfrqsPfiltersUSZ, CISTDphasePfiltersUSZ_degree,'og','markerfacec','g')
hold off
grid on
ylabel(sprintf('CI at %i%s\non the phase - degree',probaIC*100,'%'), ...
    'fontname','times','fontsize',12)
hl=legend('theoretical','emprirically');
set(hl,'fontname','times','fontsize',12)
%==================================
%
% plot(allT.TUUonUR,statURonRR.pdf,'.-')
% hold on
% plot(allT.TUUonUR,statUUonUR.pdf,'.-r')
% hold off

HorizontalSize = 16;
VerticalSize   = 20;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
    set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');

printdirectory  = ' ../../figures/';

fileprint = sprintf('%sSTDandCIonSensor%iyear%smonth%sday%s.eps',...
    printdirectory,ihc,filenameonly(7:10),filenameonly(16:17),...
    filenameonly(21:22));

fileprintepscmd = sprintf('print -depsc -loose %s',fileprint);
fileeps2pdfcmd  = sprintf('!epstopdf %s',fileprint);
filermcmd       = sprintf('!rm %s',fileprint);

%
%     eval(fileprintepscmd)
%     eval(fileeps2pdfcmd)
%     eval(filermcmd)
