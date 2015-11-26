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
directoryinputresults = '../AAresultswithFB98ter/';

%======== to compute theoretical integrals
allT.TUUonUR    = linspace(0.6,2,100);
allT.TURonRR    = linspace(0.6,2,100);
allT.MSC        = linspace(0.6,1,100);
allT.phase_rd   = linspace(-pi,pi,100);
% critical region at +/-1 sigma
alphaCI         = 0.1;
randdrawnumber  = 10;

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
    alltheoreticalmodstd         = theoreticalmodstd(:,remainindex);
    alltheoreticalphasestd_rad   = theoreticalphasestd_rad(:,remainindex);
    
    %%
    permutenbmats                = randperm(length(remainindex));
    indrandomlychosen            = permutenbmats(1:randdrawnumber);
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
    alltheoreticalmodstd         = ...
        alltheoreticalmodstd(:,indrandomlychosen);
    alltheoreticalphasestd_rad   = ...
        alltheoreticalphasestd_rad(:,indrandomlychosen);

    STDmodRatioPfilters_ave      = nanmedian(allSTDmodRatioSupPfilters,2);
    STDphaseRatioPfilters_ave    = nanmedian(allSTDphaseRatioSupPfilters,2);
    
    %== sort in increasing order, S as sort 
    [allfrqsPfiltersS, indS]     = sort(allfrqsPfilters);
    allRatioPfiltersS            = allRatioSupPfilters(indS,:);
    allmeanMSCcstPfiltersS       = allmeanMSCcstPfilters(indS,:);
    STDmodRatioPfilters_aveS     = STDmodRatioPfilters_ave(indS);
    STDphaseRatioPfilters_aveS   = STDphaseRatioPfilters_ave(indS);
    nbofvaluesoverthresholdS     = nbofvaluesoverthreshold(indS,:);
    allScpPfiltersS              = allScpPfilters(:,indS,:);
    alltheoreticalmodstdS        = alltheoreticalmodstd(indS,:);
    alltheoreticalphasestd_radS  = alltheoreticalphasestd_rad(indS,:);
    
    %== unique frequency value, U as unique
    [allfrqsPfiltersUS, indU]    = unique(allfrqsPfiltersS);
    allRatioPfiltersUS           = allRatioPfiltersS(indU,:);
    allmeanMSCcstPfiltersUS      = allmeanMSCcstPfiltersS(indU,:);
    STDmodRatioPfilters_aveUS    = STDmodRatioPfilters_aveS(indU);
    STDphaseRatioPfilters_aveUS  = STDphaseRatioPfilters_aveS(indU);
    nbofvaluesoverthresholdUS    = nbofvaluesoverthresholdS(indU,:);
    allScpPfiltersUS             = allScpPfiltersS(:,indU,:);
    alltheoreticalmodstdUS       = alltheoreticalmodstdS(indU,:);
    alltheoreticalphasestd_radUS = alltheoreticalphasestd_radS(indU,:);

    %== without 0 frequency values
    indZ                         = find(not(allfrqsPfiltersUS==0));
    allfrqsPfiltersUSZ           = allfrqsPfiltersUS(indZ);
    RatioPfiltersUSZ             = allRatioPfiltersUS(indZ,:);
    allmeanMSCcstPfiltersUSZ     = allmeanMSCcstPfiltersUS(indZ,:);
    STDmodRatioPfilters_aveUSZ   = STDmodRatioPfilters_aveUS(indZ);
    STDphaseRatioPfilters_aveUSZ = STDphaseRatioPfilters_aveUS(indZ);
    nbofvaluesoverthresholdUSZ   = nbofvaluesoverthresholdUS(indZ,:);
    allScpPfiltersUSZ            = allScpPfiltersS(:,indZ,:);
    alltheoreticalmodstdUSZ      = alltheoreticalmodstdUS(indZ,:);
    alltheoreticalphasestdUSZ_rd = alltheoreticalphasestd_radUS(indZ,:);

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
% this part performs the theoretical values of the STDs using an averaging
% for the spectral matrix SCP_ifq. Results could be more accurate. It is
% indexed by 1: CImodtheoUSZ_1 and CIphasetheoUSZ_1_degree.
%     
%====
NaverageFFTs        = filtercharact(1).ratioDFT2SCP;
LfqsUSZ             = length(allfrqsPfiltersUSZ);
% spectral matrix averaged on the selected number of records
% given by randdrawnumber for the full band of interest
meanSPsUSZ          = squeeze(nanmean(allScpPfiltersUSZ,3));

[STDmodtheoUSZ,STDphasetheoUSZ_rd] = ...
    averageSTDcomputation(allT,NaverageFFTs, meanSPsUSZ, ...
    allfrqsPfiltersUSZ, alphaCI);
sumnbofvaluesoverthresholdUSZ = sum(nbofvaluesoverthresholdUSZ,2);

% CIPfiltersUSZ = alphaCI * STDmodPfiltersUSZ ./ sumnbofvaluesoverthresholdUSZ;
CImodPfilters_aveUSZ = alphaCI * STDmodRatioPfilters_aveUSZ ./ ...
    sumnbofvaluesoverthresholdUSZ;

CImodtheoUSZ = alphaCI * STDmodtheoUSZ ./ sumnbofvaluesoverthresholdUSZ;
% CImodtheoUSZ_2 = alphaCI * nanmean(alltheoreticalmodstdUSZ,2) ./ ...
%     sumnbofvaluesoverthresholdUSZ;


% CISTDphasePfiltersUSZ_degree       = (180/pi) * alphaCI * STDphasePfiltersUSZ_rd./ ...
%     sumnbofvaluesoverthresholdUSZ;
CIphasePfilters_aveUSZ_degree   = (180/pi) * alphaCI * STDphaseRatioPfilters_aveUSZ ./ ...
    sumnbofvaluesoverthresholdUSZ;

CIphasetheoUSZ_degree = (180/pi) * alphaCI * STDphasetheoUSZ_rd ./ sumnbofvaluesoverthresholdUSZ;
% CIphasetheoUSZ_2_degree = (180/pi) * alphaCI * nanmean(alltheoreticalphasestdUSZ_rd,2) ./ ...
%     sumnbofvaluesoverthresholdUSZ;

%%
figure(2)
subplot(311)
loglog(allfrqsPfiltersUSZ,sumnbofvaluesoverthresholdUSZ,'ob','markerfacec','b')
title(sprintf('IS26 -  sensor H%i, MSC threshold = %4.2f\nday number = %i',...
    ihc, MSCthreshold, 2*randdrawnumber),'fontname','times','fontsize',14)
ylabel('counts above the threshold','fontname','times','fontsize',12)
grid on
subplot(312)
loglog(allfrqsPfiltersUSZ, CImodtheoUSZ, 'ok','markerfacec','k')
hold on
loglog(allfrqsPfiltersUSZ, CImodPfilters_aveUSZ,'or','markerfacec','r')
hold off
grid on
ylabel(sprintf('CI at %i%s\non the ratio module',(1-alphaCI)*100,'%'), ...
    'fontname','times','fontsize',12)
hl=legend('theoretical','emprirically');
set(hl,'fontname','times','fontsize',12)
posl = get(hl,'position');
set(hl,'position',[0.2 posl(2:4)])

subplot(313)
loglog(allfrqsPfiltersUSZ, CIphasetheoUSZ_degree,'ok','markerfacec','k')
hold on
loglog(allfrqsPfiltersUSZ, CIphasePfilters_aveUSZ_degree,'or','markerfacec','r')
% hold on
% loglog(allfrqsPfiltersUSZ, CISTDphasePfiltersUSZ_degree,'og','markerfacec','g')
hold off
grid on
ylabel(sprintf('CI at %i%s\non the ratio phase - degree',(1-alphaCI)*100,'%'), ...
    'fontname','times','fontsize',12)
hl=legend('theoretical','emprirically');
set(hl,'fontname','times','fontsize',12)
posl = get(hl,'position');
set(hl,'position',[0.2 posl(2:4)])
%==================================
%
% plot(allT.TUUonUR,statURonRR.pdf,'.-')
% hold on
% plot(allT.TUUonUR,statUUonUR.pdf,'.-r')
% hold off

HorizontalSize = 16;
VerticalSize   = 22;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
    set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');

printdirectory  = ' ../../figures/';

fileprint = sprintf('%sSTDandCIonSensor%iyear%smonth%sday%snumber%i.eps',...
    printdirectory,ihc,filenameonly(7:10),filenameonly(16:17),...
    filenameonly(21:22),randdrawnumber);

fileprintepscmd = sprintf('print -depsc -loose %s',fileprint);
fileeps2pdfcmd  = sprintf('!epstopdf %s',fileprint);
filermcmd       = sprintf('!rm %s',fileprint);

%
    eval(fileprintepscmd)
    eval(fileeps2pdfcmd)
    eval(filermcmd)
