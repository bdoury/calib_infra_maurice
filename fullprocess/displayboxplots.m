%========================== displaySensorRatio.m =========================
% this program reads the ratios SUT/SREF estimated by spectral approach
% and stored in the drectory AAresults. That consists on 8 files
% for the 8 sensors of IS26. Each file consists the estimate ratios
% on several days of records. The
%
% this program calls the program theoreticalSUTresponse4IS26.m
%
clear

addpath ZZtoolbox/

sensor_UT = 'I26DE_BDF_RSP_2015134_MB3';

for ihc = 1:5
    
    switch ihc
        case 1
            coeffsens=1.037;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 2
            coeffsens=1.065;%.1;
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
            coeffsens=1.0;
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
    
    %     %======= fit the curve
    %     segmentsnumber = 2;
    %     polydegrees    = (0:1:12);
    %     logtrain.flag  = 1;
    %     logtrain.N     = 700;
    %     logfit.flag    = 1;
    %     logfit.N       = 200;
    %
    %     [FreqFitabs_Hz,absRfit] = smoothpolyLL(allfrqsPfilters,...
    %         absallRatioPfilters_ave, ...
    %         segmentsnumber,polydegrees,logtrain,logfit);
    %
    %     [FreqFitphase_Hz,phaseRfit_deg] = smoothpolyLL(allfrqsPfilters,...
    %         angleallRatioPfilters_ave_deg,...
    %         segmentsnumber,polydegrees,logtrain,logfit);
    
    %=======
    
    [allfrqsPfiltersU, inda] = unique(allfrqsPfilters);
    allRatioPfiltersU        = allRatioPfilters(inda,:);
    allRatioPfiltersU        = allRatioPfiltersU(not(allfrqsPfiltersU==0),:);
    allfrqsPfiltersU         = allfrqsPfiltersU(not(allfrqsPfiltersU==0));
 
    
    
    [p_total_NRS_sensor, freq_vector, p_total_NRS, TF_ref_sensor, TFsensor4freqRatio] = ...
        theoreticalSUTresponse4IS26(allfrqsPfiltersU, sensor_UT, ref_sensor, 'nofir');
    
    absestim = coeffsens * abs(allRatioPfiltersU) .* ...
        (abs(TFsensor4freqRatio)*ones(1,size(allRatioPfiltersU,2)));
    
    figure(ihc)
    
    %================================================
    subplot(121)
    boxplot(absestim','position',allfrqsPfiltersU,'symbol','','whisker',0);
    set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
        'xticklabel',[0.001 0.01 0.1 1 10])
    set(gca,'fontname','times','fontsize',10)
    set(gca,'xlim',[0.008 6])
    set(gca,'ylim',[0.8 1.2])
    grid on
    xlabel('frequency - Hz')
    ylabel('estimated gain')
        
    %=============== dipslay
    hold on
    semilogx(freq_vector, abs(p_total_NRS_sensor), 'r');
    semilogx(freq_vector, abs(p_total_NRS_sensor)*0.95, 'r--','linew',1.5);
    semilogx(freq_vector, abs(p_total_NRS_sensor)*1.05, 'r--','linew',1.5);
    hold off
    title(sprintf('IS26 -  sensor #%i, threshold = %4.2f\nday number = %i',...
        ihc, MSCthreshold, 2*doubledaynumber),'fontname','times','fontsize',12)

    %========================== PHASE =========
    anglestime_deg = angle(allRatioPfiltersU)*180/pi;
    if ihc<=4
        anglestime_deg=-anglestime_deg;
    end
    
    anglestime_deg = anglestime_deg + ...
        (angle(TFsensor4freqRatio)*ones(1,size(anglestime_deg,2)))*180/pi;
    subplot(122)
    boxplot(anglestime_deg','position',allfrqsPfiltersU,'symbol','','whisker',0);
    set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
        'xticklabel',[0.001 0.01 0.1 1 10])
    set(gca,'fontname','times','fontsize',10)
    set(gca,'xlim',[0.008 6])
    set(gca,'ylim',[-10 10])
    grid on
    xlabel('frequency - Hz')
    ylabel('estimated phase - degree')
    
    
    %=============== dipslay
    hold on
    semilogx(freq_vector, angle(p_total_NRS_sensor)*180/pi, 'r');
    semilogx(freq_vector, (angle(p_total_NRS_sensor)*180/pi)-5, 'r--','linew',1.5);
    semilogx(freq_vector, (angle(p_total_NRS_sensor)*180/pi)+5, 'r--','linew',1.5);
    hold off

    
    %==============================================================
    HorizontalSize = 22;
    VerticalSize   = 14;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a3');
    set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    set(gcf,'color', [1,1,0.92]);
    set(gcf, 'InvertHardCopy', 'off');
    
    printdirectory  = ' ../slidesITW2015/';
    fileprintepscmd = sprintf('print -depsc -loose %s3monthsonIS26SUTboxplot%i.eps',printdirectory,ihc);
    fileeps2pdfcmd  = sprintf('!epstopdf %s3monthsonIS26SUT%i.eps',printdirectory,ihc);
    filermcmd       = sprintf('!rm %s3monthsonIS26SUT%i.eps',printdirectory,ihc);
    %
      eval(fileprintepscmd)
    %     eval(fileeps2pdfcmd)
    %     eval(filermcmd)
    
    %     fileprintpngcmd = sprintf('print -dpng -loose ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
    %     fileprintpngcmd = sprintf('print -dpdf -loose threemonthsonIS26SUT%i',ihc);
    %     eval(fileprintpngcmd)
    
    
    
    
    
end
