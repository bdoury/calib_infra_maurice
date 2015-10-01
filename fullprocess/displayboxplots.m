%========================== displaySensorRatio.m =========================
% this program reads the ratios SUT/SREF estimated by spectral approach
% and stored in a drectory as AAresults. That consists on 8 files
% for the 8 sensors of IS26. Each file consists the estimate ratios
% on several days of records. 
%
% this program calls the program theoreticalSUTresponse4IS26.m
%
clear
addpath ZZtoolbox/
addpath ZZtoolbox/00gabrielson
sensor_UT = 'I26DE_BDF_RSP_2015134_MB3';
for ihc = 1:5
    switch ihc
        case 1
            coeffsens=1.048;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 2
            coeffsens=1.08;%.1;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 3
            coeffsens=1.065;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 4
            coeffsens=1.07;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 5
            coeffsens=1.08;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 6
            coeffsens=0.97;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB3';
        case 7
            coeffsens=0.97;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB3';
        case 8
            coeffsens=0.99;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB3';
    end
    
    comload = sprintf('load AAresultswithFB/resultssta26sensor%i',ihc);
    eval(comload);
    Dstart = str2double(fileswithdotmat(1).name(13:15));
    Dend   = str2double(fileswithdotmat(length(fileswithdotmat)).name(13:15));
    if 1
        doubledaynumber = (Dend-Dstart+3)/2;
    else
        doubledaynumber = 10;
        permutenbmats = randperm(nbmats);
        allRatioPfilters = allRatioPfilters(:,permutenbmats(1:doubledaynumber));
    end
    
    allRatioPfilters_ave          = nanmean(allRatioPfilters,2);
    absallRatioPfilters_ave       = abs(allRatioPfilters_ave);
    angleallRatioPfilters_ave_deg = angle(allRatioPfilters_ave)*180/pi;
    
    %=======
    if 1
        %======= fit the curve
        segmentsnumber = 2;
        polydegrees    = (0:1:8);
        logtrain.flag  = 1;
        logtrain.N     = 115;
        logfit.flag    = 1;
        logfit.N       = 115;
    else
        %======= fit the curve
        segmentsnumber = 2;
        polydegrees    = (0:1:8);
        logtrain.flag  = 1;
        logtrain.N     = 50;
        logfit.flag    = 1;
        logfit.N       = 200;
    end
    reducefreqrange = (1:fix(0.97*length(allfrqsPfilters)));

    
    [allfrqsPfiltersU, inda] = unique(allfrqsPfilters);
    allRatioPfiltersU        = allRatioPfilters(inda,:);
    allmeanMSCcstPfiltersU   = allmeanMSCcstPfilters(inda,:);
    nbofvaluesoverthresholdU = nbofvaluesoverthreshold(inda,:);
    
    allRatioPfiltersUZ        = allRatioPfiltersU(not(allfrqsPfiltersU==0),:);
    allmeanMSCcstPfiltersUZ   = allmeanMSCcstPfiltersU(not(allfrqsPfiltersU==0),:);
    allfrqsPfiltersUZ         = allfrqsPfiltersU(not(allfrqsPfiltersU==0));
    nbofvaluesoverthresholdUZ = nbofvaluesoverthresholdU(not(allfrqsPfiltersU==0),:);
    
    meanallRatioPfiltersUZ    = nanmean(allRatioPfiltersUZ,2);
    
    stdallRatioPfiltersUZ    = nanstd(allRatioPfiltersUZ,[],2);
    
    ICallRatioPfiltersUZ     = stdallRatioPfiltersUZ ./ sqrt(sum(nbofvaluesoverthresholdUZ,2));
    
    
    N_freq_vector = 400;
    freq_vector = logspace(log10(0.001),log10(30),N_freq_vector) .';
    [p_total_NRS_sensor, p_total_NRS, TF_ref_sensor, TFsensor4freqRatio] = ...
        HCP_acoustical(freq_vector, allfrqsPfiltersUZ, sensor_UT, ref_sensor, 'nofir');
    
    absestim = coeffsens * abs(meanallRatioPfiltersUZ) .* abs(TFsensor4freqRatio);
    
    
%     [FreqFitabs_Hz,absRfit] = smoothpolyLL(allfrqsPfiltersUZ, ...
%         coeffsens * abs(meanallRatioPfiltersUZ),...
%         segmentsnumber,polydegrees,logtrain,logfit);
    
    figure(ihc)
    clf
    %================================================
     subplot(211)
    semilogx(allfrqsPfiltersUZ,20*log10(absestim),'.-k'),
%     hold on
%         semilogx(ones(2,1)*allfrqsPfiltersUZ',...
%             [absestim'-ICallRatioPfiltersUZ';...
%             absestim'+ICallRatioPfiltersUZ'],'.-b')

%     boxplot(absestim','position',allfrqsPfiltersUZ,'symbol','','whisker',0);
    set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
        'xticklabel',[0.001 0.01 0.1 1 10])
    set(gca,'fontname','times','fontsize',14)
    set(gca,'xlim',[0.008 10])
    set(gca,'ylim',[-5 5])
    grid on
%     xlabel('frequency [Hz]')
    ylabel('Amplitude [dB]','fontname','times','fontsize',14)
    hold on
    plot([1 1]*0.02,[-45 45],'m','linew',2)
    plot([1 1]*4,[-45 45],'m','linew',2)
    hold off

    %=============== dipslay
    hold on
%      semilogx(freq_vector, 20*log10(abs(p_total_NRS_sensor)), 'r');
    semilogx(freq_vector, 20*log10(abs(p_total_NRS_sensor)*0.95), 'r--','linew',1.5);
    semilogx(freq_vector, 20*log10(abs(p_total_NRS_sensor)*1.05), 'r--','linew',1.5);
    hold off
%     title(sprintf('IS26 -  sensor H%i, MSC threshold = %4.2f\nday number = %i',...
%         ihc, MSCthreshold, 2*doubledaynumber),'fontname','times','fontsize',14)
    title(sprintf('IS26 -  sensor H%i\ndashed line: +/-5%s for amplitude, +/- 5 degrees for phase', ihc,'%'),...
        'fontname','times','fontsize',14)

    %========================== PHASE =========
    anglestime_deg = -angle(meanallRatioPfiltersUZ)*180/pi;
%     if ihc<=4
%         anglestime_deg=-anglestime_deg;
%     end
    
%     anglestime_deg = anglestime_deg + angle(TFsensor4freqRatio)*180/pi;
    subplot(212)
    semilogx(allfrqsPfiltersUZ,unwrap(anglestime_deg*pi/180)*180/pi,'.-k'),

%     boxplot(anglestime_deg','position',allfrqsPfiltersUZ,'symbol','','whisker',0);
    set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
        'xticklabel',[0.001 0.01 0.1 1 10])
    set(gca,'fontname','times','fontsize',14)
    set(gca,'xlim',[0.008 10])
    set(gca,'ylim',[-20 20])
    grid on
    xlabel('frequency [Hz]')
    ylabel('Phase [deg]')
    
    hold on
    plot([1 1]*0.02,[-45 45],'m','linew',2)
    plot([1 1]*4,[-45 45],'m','linew',2)
    hold off
    %=============== dipslay
    hold on
%     semilogx(freq_vector, -unwrap(angle(p_total_NRS))*180/pi, 'r');
    semilogx(freq_vector, -unwrap(angle(p_total_NRS))*180/pi-5, 'r--','linew',1.5);
    semilogx(freq_vector, -unwrap(angle(p_total_NRS))*180/pi+5, 'r--','linew',1.5);
    hold off
ht = text(0.08, -16,'IMS passband');
set(ht,'color','m','fontsize',14,'fontname','times')
    
    %==============================================================
    HorizontalSize = 16;
    VerticalSize   = 10;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a3');
    set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    set(gcf,'color', [1,1,0.92]);
    set(gcf, 'InvertHardCopy', 'off');
    
    printdirectory  = ' ../slidesITW2015/';
    fileprintepscmd = sprintf('print -depsc -loose %s3monthsonIS26SUTboxplot%i.eps',printdirectory,ihc);
    fileeps2pdfcmd  = sprintf('!epstopdf %s3monthsonIS26SUTboxplot%i.eps',printdirectory,ihc);
    filermcmd       = sprintf('!rm %s3monthsonIS26SUTboxplot%i.eps',printdirectory,ihc);
    %
%       eval(fileprintepscmd)
%         eval(fileeps2pdfcmd)
%         eval(filermcmd)
    
    %     fileprintpngcmd = sprintf('print -dpng -loose ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
    %     fileprintpngcmd = sprintf('print -dpdf -loose threemonthsonIS26SUT%i',ihc);
    %     eval(fileprintpngcmd)
    
    
    
    
    
end
