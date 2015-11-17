%========================== displaySUTresponse.m =========================
% this program reads the ratios SUT/SREF estimated by spectral approach
% and stored in a drectory as AAresults. That consists on 8 files
% for the 8 sensors of IS26. Each file consists the estimate ratios
% on several days of records.
%=========
% this program calls the function:
%         - HCP_acoustical.m
% using both ZZtoolbox/00gabrielson and ../ZZtoolbox/00benoit
%
%=========================================================================
clear
addpath ../ZZtoolbox/
addpath ../ZZtoolbox/00gabrielson
addpath ../ZZtoolbox/00benoit

%==== this directory contains the parameters evalauted by the
% program estimationwithFB.m
directoryinputresults = '../AAresultswithFB98/';
% directoryinputresults = '../AAresultswithFBter/';
% directoryinputresults = '../AAresultswithFBbis/';

sensor_UT = 'I26DE_BDF_RSP_2015134_MB3';
saveflag = 0;
    remainindex = [20:40];

for ihc = 2
    comload = sprintf('load %sresultssta26sensor%i',directoryinputresults,ihc);
    numfig = ihc+100;
    switch ihc
        case 1
            coeffsens=1.04;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 2
            coeffsens=1.1;%.1;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 3
            coeffsens=1.065;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 4
            coeffsens=1.07;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 5
            coeffsens=1.1;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 6
            coeffsens=1;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB3';
        case 7
            coeffsens=0.97;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB3';
        case 8
            coeffsens=1.02;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB3';
    end
    eval(comload);
    %%
    Dstart = str2double(fileswithdotmat(1).name(13:15));
    Dend   = str2double(fileswithdotmat(length(fileswithdotmat)).name(13:15));
    if 1
        doubledaynumber = (Dend-Dstart+3)/2;
    else
        doubledaynumber = 10;
        permutenbmats = randperm(nbmats);
        allRatioPfilters = allRatioSupPfilters(:,permutenbmats(1:doubledaynumber));
    end
    
    STDmodRatioPfilters_ave    = nanmean(allSTDmodRatioSupPfilters,2);
    STDphaseRatioPfilters_ave  = nanmean(allSTDphaseRatioSupPfilters,2);
    
    %== sort in increasing order
    [allfrqsPfiltersS, inds]       = sort(allfrqsPfilters);
    STDmodRatioPfilters_aveS       = STDmodRatioPfilters_ave(inds);
    STDphaseRatioPfilters_aveS     = STDphaseRatioPfilters_ave(inds);
    allRatioPfiltersS              = allRatioSupPfilters(inds,remainindex);
    allmeanMSCcstPfiltersS         = allmeanMSCcstPfilters(inds,remainindex);
    nbofvaluesoverthresholdS       = nbofvaluesoverthreshold(inds,remainindex);

    %== unique
    [allfrqsPfiltersUS, inda]      = unique(allfrqsPfiltersS);
    STDmodRatioPfilters_aveUS      = STDmodRatioPfilters_aveS(inda);
    STDphaseRatioPfilters_aveUS    = STDphaseRatioPfilters_aveS(inda);
    
    allRatioPfiltersUS             = allRatioPfiltersS(inda,:);
    allmeanMSCcstPfiltersUS        = allmeanMSCcstPfiltersS(inda,:);
    nbofvaluesoverthresholdUS      = nbofvaluesoverthresholdS(inda,:);    
    
    %== without 0
    RatioPfiltersUSZ               = allRatioPfiltersUS(not(allfrqsPfiltersUS==0),:);
    allmeanMSCcstPfiltersUSZ       = allmeanMSCcstPfiltersUS(not(allfrqsPfiltersUS==0),:);
    allfrqsPfiltersUSZ             = allfrqsPfiltersUS(not(allfrqsPfiltersUS==0));
    nbofvaluesoverthresholdUSZ     = nbofvaluesoverthresholdUS(not(allfrqsPfiltersUS==0),:);
    STDmodRatioPfilters_aveUSZ     = STDmodRatioPfilters_aveUS(not(allfrqsPfiltersUS==0));
    STDphaseRatioPfilters_aveUSZ   = STDphaseRatioPfilters_aveUS(not(allfrqsPfiltersUS==0));
    
    %======
    modRatioPfiltersUSZ            = abs(RatioPfiltersUSZ);
    meanmodRatioPfiltersUSZ        = nanmean(modRatioPfiltersUSZ,2);    
    STDmodPfiltersUSZ              = nanstd(modRatioPfiltersUSZ,[],2);
    
    phaseRatioPfiltersUSZ_rd       = angle(RatioPfiltersUSZ);
    meanphasePfiltersUSZ_rd        = nanmean(phaseRatioPfiltersUSZ_rd,2);    
    STDphasePfiltersUSZ_rd         = nanstd(phaseRatioPfiltersUSZ_rd,[],2);

    ICallRatioPfiltersUSZ          = STDmodPfiltersUSZ ./ ...
        sqrt(sum(nbofvaluesoverthresholdUSZ,2));
    
    ICallRatioPfiltersUSZbis       = STDmodRatioPfilters_aveUSZ ./ ...
        sqrt(sum(nbofvaluesoverthresholdUSZ,2));
    
    N_freq_vector = 300;
    freq_vector   = logspace(log10(0.001),log10(30),N_freq_vector) .';
    [p_total_NRS_sensor, p_total_NRS, TF_ref_sensor, TFsensor4freqRatio] = ...
        HCP_acoustical(freq_vector, allfrqsPfiltersUSZ, sensor_UT, ref_sensor, 'nofir');
        
    %================================
    absestimwithcorrect = coeffsens * meanmodRatioPfiltersUSZ .* abs(TFsensor4freqRatio);   
    figure(numfig)
    clf
    %================================================
    subplot(211)
    semilogx(allfrqsPfiltersUSZ,20*log10(absestimwithcorrect),'k','linew',1.5),
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
    %     boxplot(absestim','position',allfrqsPfiltersUZ,'symbol','','whisker',0);
    %         set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
    %             'xticklabel',[0.001 0.01 0.1 1 10])
    
    %         set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
    %             'xticklabel',[0.001 0.01 0.1 1 10])
    grid on
    ylabel('Amplitude [dB]','fontname','times','fontsize',14)
    hold on
    plot([1 1]*0.0205,[-45 45],'r','linew',2)
    plot([1 1]*3.98,[-45 45],'r','linew',2)
    hold off
    
    %=============== dipslay
    hold on
    %      semilogx(freq_vector, 20*log10(abs(p_total_NRS_sensor)), 'r');
    semilogx(freq_vector, 20*log10(abs(p_total_NRS_sensor)*0.95), 'r--','linew',1.5);
    semilogx(freq_vector, 20*log10(abs(p_total_NRS_sensor)*1.05), 'r--','linew',1.5);
    hold off
    title(sprintf('IS26 -  sensor H%i, MSC threshold = %4.2f\nday number = %i',...
        ihc, MSCthreshold, 2*doubledaynumber),'fontname','times','fontsize',14)
    %         title(sprintf('IS26 -  sensor H%i\ndashed line: +/-5%s for amplitude, +/- 5 degrees for phase', ihc,'%'),...
    %             'fontname','times','fontsize',14)
    
    set(gca,'fontname','times','fontsize',14)
    set(gca,'xlim',[0.01 5])
    set(gca,'ylim',4*[-1 1])
    set(gca,'xtickLabel',[])
    %              xlabel('frequency [Hz]')
    
    set(gca,'position',[0.1300    0.5056    0.7750    0.3559])
    
    ht = text(0.022,-3.4,'0.02 Hz');
    set(ht,'color','r','fontsize',14,'fontname','times')
    ht = text(2.1,-3.4,'4 Hz');
    set(ht,'color','r','fontsize',14,'fontname','times')
    ht = text(0.14, -3.4,'IMS passband');
    set(ht,'color','r','fontsize',14,'fontname','times')
    
    %========================== PHASE =========
    
    anglewithcorrect_rd = meanphasePfiltersUSZ_rd + angle(TFsensor4freqRatio);
    angltheo_rd   = angle(p_total_NRS)+ angle(TF_ref_sensor);
    
    subplot(212)
    
    %         semilogx(allfrqsPfiltersUZ,unwrap(anglestime_rd)*180/pi,'.-k'),
    
    %     boxplot(anglestime_deg','position',allfrqsPfiltersUZ,'symbol','','whisker',0);
    %         set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
    %             'xticklabel',[0.001 0.01 0.1 1 10])
    semilogx(allfrqsPfiltersUSZ,unwrap(anglewithcorrect_rd)*180/pi,'k','linew',1.5),
    
    %     boxplot(anglestime_deg','position',allfrqsPfiltersUZ,'symbol','','whisker',0);
    %         set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
    %             'xticklabel',[0.001 0.01 0.1 1 10])
    set(gca,'fontname','times','fontsize',14)
    
    grid on
    xlabel('frequency [Hz]')
    ylabel('Phase [deg]')
    
    hold on
    
    plot([1 1]*0.0205,[-45 45],'r','linew',2)
    plot([1 1]*3.98,[-45 45],'r','linew',2)
    hold off
    %=============== dipslay
    hold on
    %     semilogx(freq_vector, -unwrap(angltheo_rd)*180/pi, 'r');
    semilogx(freq_vector, unwrap(angltheo_rd)*180/pi-5, 'r--','linew',1.5);
    semilogx(freq_vector, unwrap(angltheo_rd)*180/pi+5, 'r--','linew',1.5);
    hold off
    
    set(gca,'fontname','times','fontsize',14)
    set(gca,'xlim',[0.01 5])
    set(gca,'ylim',40*[-1 1])
    xlabel('frequency [Hz]')
    
    set(gca,'position',[ 0.1300    0.1328    0.7750    0.3333])
    %==============================================================
    HorizontalSize = 16;
    VerticalSize   = 10;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a3');
    set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    %         set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    set(gcf,'color', [1,1,0.92]);
    set(gcf, 'InvertHardCopy', 'off');
    
    printdirectory  = ' ../slidesITW2015/';
    fileprintepscmd = sprintf('print -depsc -loose %s3monthsonIS26SUTboxplot%i.eps',printdirectory,ihc);
    fileeps2pdfcmd  = sprintf('!epstopdf %s3monthsonIS26SUTboxplot%i.eps',printdirectory,ihc);
    filermcmd       = sprintf('!rm %s3monthsonIS26SUTboxplot%i.eps',printdirectory,ihc);
    
end

if saveflag
    eval(fileprintepscmd)
    eval(fileeps2pdfcmd)
    eval(filermcmd)
end
%%
% 
% semilogx(allfrqsPfilters, sum(nbofvaluesoverthreshold,2),'.')
% hold on
% for ip=1:Pfilter
%     plot(allfrqsFFT_Hz{ip}(idipinf(ip))*ones(2,1),[0 20000],allcolors(ip,1))
%     plot(allfrqsFFT_Hz{ip}(idipsup(ip))*ones(2,1),[0 20000],allcolors(ip,1))
%     pause
% end
% hold off
%%
% NaverageFFTs = filtercharact(1).ratioDFT2SCP;
% allT.TUUonUR = linspace(0.6,3,100);
% allT.TURonRR = linspace(0.6,3,100);
% allT.MSC     = linspace(0.5,1,100);
% allT.phase   = linspace(0,2*pi,100);
%
% [statUUonUR, statURonRR, statMSC]    = ...
%     statsRatiosHbis(allT, HH, NaverageFFTs, 0.7);

