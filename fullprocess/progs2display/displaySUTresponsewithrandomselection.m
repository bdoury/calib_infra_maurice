%======== displaySUTresponsewithrandomselection.m =========================
% this program reads the ratios SUT/SREF estimated by spectral approach
% and stored in a drectory as AAresults. That consists on 8 files
% for the 8 sensors of IS26. Each file consists the estimate ratios
% on several days of records.
% We do not correct by sensor/WNRS responses
% We randomly choose the files to be averaged in such a way
% we observe a certain stability and avoid outliers by trimmed means
%=========================================================================
clear
addpath ../ZZtoolbox/
directorysignals    = '../../../../AAdataI26calib/';

%==== this directory contains the parameters evalauted by the
% program estimationwithFB.m
directoryinputresults = '../AAresultswithFB98/';

saveflag                = 0;
randomlydoubledaynumber = 15;

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
            remainindex = [1:61 62 63:70]; % 2015/10/13
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
    
    
    permutenbmats     = randperm(length(remainindex));
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
    
    %
    %     STDmodRatioPfilters_ave        = nanmedian(allSTDmodRatioSupPfilters,2);
    %     STDphaseRatioPfilters_ave      = nanmedian(allSTDphaseRatioSupPfilters,2);
    
    %== sort in increasing order
    [allfrqsPfiltersS, inds]       = sort(allfrqsPfilters);
    allRatioPfiltersS              = allRatioSupPfilters(inds,:);
    allmeanMSCcstPfiltersS         = allmeanMSCcstPfilters(inds,:);
    %     STDmodRatioPfilters_aveS       = STDmodRatioPfilters_ave(inds);
    %     STDphaseRatioPfilters_aveS     = STDphaseRatioPfilters_ave(inds);
    nbofvaluesoverthresholdS       = nbofvaluesoverthreshold(inds,:);
    
    %== unique frequency value
    [allfrqsPfiltersUS, inda]      = unique(allfrqsPfiltersS);
    allRatioPfiltersUS             = allRatioPfiltersS(inda,:);
    allmeanMSCcstPfiltersUS        = allmeanMSCcstPfiltersS(inda,:);
    %     STDmodRatioPfilters_aveUS      = STDmodRatioPfilters_aveS(inda);
    %     STDphaseRatioPfilters_aveUS    = STDphaseRatioPfilters_aveS(inda);
    nbofvaluesoverthresholdUS      = nbofvaluesoverthresholdS(inda,:);
    
    %== without 0 frequency values
    allfrqsPfiltersUSZ             = allfrqsPfiltersUS(not(allfrqsPfiltersUS==0));
    RatioPfiltersUSZ               = allRatioPfiltersUS(not(allfrqsPfiltersUS==0),:);
    allmeanMSCcstPfiltersUSZ       = allmeanMSCcstPfiltersUS(not(allfrqsPfiltersUS==0),:);
    %     STDmodRatioPfilters_aveUSZ     = STDmodRatioPfilters_aveUS(not(allfrqsPfiltersUS==0));
    %     STDphaseRatioPfilters_aveUSZ   = STDphaseRatioPfilters_aveUS(not(allfrqsPfiltersUS==0));
    nbofvaluesoverthresholdUSZ     = nbofvaluesoverthresholdUS(not(allfrqsPfiltersUS==0),:);
    
    %====== absolute and arg of the ratios
    modRatioPfiltersUSZ            = abs(RatioPfiltersUSZ);
    phaseRatioPfiltersUSZ_rd       = angle(RatioPfiltersUSZ);
    %====== averaging by TRIMMEAN to avoid outliers
    trimmeanmodRatioPfiltersUSZ    = trimmean(modRatioPfiltersUSZ,30,2);
    meanmodRatioPfiltersUSZ        = nanmean(modRatioPfiltersUSZ,2);
    meanphasePfiltersUSZ_rd        = trimmean(phaseRatioPfiltersUSZ_rd,30,2);
    %     %====== STDs
    %     STDmodPfiltersUSZ              = nanstd(modRatioPfiltersUSZ,[],2);
    %     STDphasePfiltersUSZ_rd         = nanstd(phaseRatioPfiltersUSZ_rd,[],2);
    %
    %     ICallRatioPfiltersUSZ          = STDmodPfiltersUSZ ./ ...
    %         sqrt(sum(nbofvaluesoverthresholdUSZ,2));
    %
    %     ICallRatioPfiltersUSZbis       = STDmodRatioPfilters_aveUSZ ./ ...
    %         sqrt(sum(nbofvaluesoverthresholdUSZ,2));
    %
    %================================================
    %     figure(numfig)
    clf
    %     subplot(211)
    semilogx(allfrqsPfiltersUSZ,(modRatioPfiltersUSZ),'o',...
        'color',0.7*[1 1 1]);
    hold on
    semilogx(allfrqsPfiltersUSZ,(meanmodRatioPfiltersUSZ),'ko',...
        'markersize',6,'markerfaceco','r'),
    semilogx(allfrqsPfiltersUSZ,(trimmeanmodRatioPfiltersUSZ),'ko',...
        'markersize',6,'markerfaceco','b'),
    hold off
    
    %         hold on
    %             semilogx(allfrqsPfiltersUZ,...
    %                 abs(meanallRatioPfiltersUZ)-allSTDmodRatioPfilters_aveU,'.-b', ...
    %                 allfrqsPfiltersUZ,...
    %                 abs(meanallRatioPfiltersUZ)+allSTDmodRatioPfilters_aveU,'.-r')
    %             hold on
    %             semilogx(allfrqsPfiltersUZ,...
    %                 abs(meanallRatioPfiltersUZ)-ICallRatioPfiltersUZ,'.g', ...
    %                 allfrqsPfiltersUZ,...
    %                 abs(meanallRatioPfiltersUZ)+ICallRatioPfiltersUZ,'.m')
    %             hold off
    % %
    %         boxplot(modRatioPfiltersUSZ','position',allfrqsPfiltersUSZ,'symbol','','whisker',0);
    %             set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
    %                 'xticklabel',[0.001 0.01 0.1 1 10])
    
    %         set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
    %             'xticklabel',[0.001 0.01 0.1 1 10])
    grid on
    ylabel('Amplitude [dB]','fontname','times','fontsize',14)
    
    %=============== dipslay
    hold off
    title(sprintf('IS26 -  sensor H%i, MSC threshold = %4.2f\nday number = %i',...
        ihc, MSCthreshold, 2*randomlydoubledaynumber),'fontname','times','fontsize',14)
    %         title(sprintf('IS26 -  sensor H%i\ndashed line: +/-5%s for amplitude, +/- 5 degrees for phase', ihc,'%'),...
    %             'fontname','times','fontsize',14)
    
    set(gca,'fontname','times','fontsize',14)
    set(gca,'xlim',[0.01 5])
    set(gca,'ylim',[0.8 1.2])
    set(gca,'xtickLabel',[])
    %              xlabel('frequency [Hz]')
    
    
    
    %========================== PHASE =========
    if 0
        anglewithoutcorrect_rd = meanphasePfiltersUSZ_rd;
        %     set(gca,'position',[ 0.1300    0.1328    0.7750    0.3333])
        %     set(gca,'position',[0.1300    0.5056    0.7750    0.3559])
        subplot(212)
        
        %         semilogx(allfrqsPfiltersUZ,unwrap(anglestime_rd)*180/pi,'.-k'),
        
        %     boxplot(anglestime_deg','position',allfrqsPfiltersUZ,'symbol','','whisker',0);
        %         set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
        %             'xticklabel',[0.001 0.01 0.1 1 10])
        semilogx(allfrqsPfiltersUSZ,unwrap(anglewithoutcorrect_rd)*180/pi,...
            'ko','markersize',6,'markerfaceco','R')
        
        %     boxplot(anglestime_deg','position',allfrqsPfiltersUZ,'symbol','','whisker',0);
        %         set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
        %             'xticklabel',[0.001 0.01 0.1 1 10])
        set(gca,'fontname','times','fontsize',14)
        
        grid on
        xlabel('frequency [Hz]')
        ylabel('Phase [deg]')
        %=============== dipslay
        
        set(gca,'fontname','times','fontsize',14)
        set(gca,'xlim',[0.01 5])
        set(gca,'ylim',40*[-1 1])
    end
    xlabel('frequency [Hz]')
    
    %==============================================================
    HorizontalSize = 16;
    VerticalSize   = 10;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a3');
    %         set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    set(gcf,'color', [1,1,0.92]);
    set(gcf, 'InvertHardCopy', 'off');
    
    printdirectory  = ' ../../figures/';
    
    fileprint = sprintf('%swithtrimmeanonstation%i.eps',printdirectory,ihc);
    
    fileprintepscmd = sprintf('print -depsc -loose %s',fileprint);
    fileeps2pdfcmd  = sprintf('!epstopdf %s',fileprint);
    filermcmd       = sprintf('!rm %s',fileprint);
    
end

if saveflag
    eval(fileprintepscmd)
    eval(fileeps2pdfcmd)
    eval(filermcmd)
end
