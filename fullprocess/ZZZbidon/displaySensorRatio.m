%========================== displaySensorRatio.m =========================
% this program reads the ratios SUT/SREF estimated by spectral approach
% and stored in the drectory AAresults. That consists on 8 files 
% for the 8 sensors of IS26. Each file consists the estimate ratios
% on several days of records. The 
%
%
%
clear

addpath ZZtoolbox/

for ihc = 5:5
    figure(ihc)
    comload = sprintf('load AAresults/resultssta26sensor%i',ihc);
    eval(comload);
    Dstart = str2double(filesmat(1).name(13:15));
    Dend   = str2double(filesmat(length(filesmat)).name(13:15));
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

    %======= fit the curve
    segmentsnumber = 2;
    polydegrees    = (0:1:12);
    logtrain.flag  = 1;
    logtrain.N     = 700;
    logfit.flag    = 1;
    logfit.N       = 200;

    [FreqFitabs_Hz,absRfit] = smoothpolyLL(allfrqsPfilters,...
        absallRatioPfilters_ave, ...
        segmentsnumber,polydegrees,logtrain,logfit);

    [FreqFitphase_Hz,phaseRfit_deg] = smoothpolyLL(allfrqsPfilters,...
        angleallRatioPfilters_ave_deg,...
        segmentsnumber,polydegrees,logtrain,logfit);
    
    %=======
    subplot(2,1,1)
    semilogx(allfrqsPfilters,20*log10(abs(allRatioPfilters)),'.', ...
        'color',0.9*ones(3,1))
    hold on
    semilogx(allfrqsPfilters,20*log10(absallRatioPfilters_ave),'.')
    hold on
    semilogx(FreqFitabs_Hz,20*log10(absRfit),'r')
    hold off
    set(gca,'fontname','times','fontsize',10)
    set(gca,'xlim',[0.002 4])
    set(gca,'ylim',[-2 2])
    grid on
    xlabel('frequency - Hz')
    ylabel('ratio magnitude - dB')
    drawnow
    %=====
    subplot(212)
        semilogx(allfrqsPfilters,angle(allRatioPfilters)*180/pi,'.', ...
        'color',0.7*ones(3,1))
    hold on
    semilogx(allfrqsPfilters,angleallRatioPfilters_ave_deg,'.')
    hold on
    semilogx(FreqFitphase_Hz,phaseRfit_deg,'r')
    hold off
    set(gca,'fontname','times','fontsize',10)
    set(gca,'xlim',[0.002 4])
    set(gca,'ylim',[-15 15])
    grid on
    xlabel('frequency - Hz')
    ylabel('ratio phase - degree')
    drawnow

    subplot(211)
    title(sprintf('IS26 -  sensor #%i, threshold = %4.2f\nday number = %i',...
        ihc, MSCthreshold, 2*doubledaynumber),'fontname','times','fontsize',12)
    HorizontalSize = 16;
    VerticalSize   = 12;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a3');
    %         set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    set(gcf,'color', [1,1,0.92]);
    set(gcf, 'InvertHardCopy', 'off');
    
    printdirectory = ' ../../../textes/6distConjointHMSC/figures/';
    fileprintepscmd = sprintf('print -depsc -loose %s3monthsonIS26SUT%i.eps',printdirectory,ihc);
    fileeps2pdfcmd = sprintf('!epstopdf %s3monthsonIS26SUT%i.eps',printdirectory,ihc);
    filermcmd = sprintf('!rm %s3monthsonIS26SUT%i.eps',printdirectory,ihc);
% 
%     eval(fileprintepscmd)
%     eval(fileeps2pdfcmd)
%     eval(filermcmd)

%     fileprintpngcmd = sprintf('print -dpng -loose ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
%     fileprintpngcmd = sprintf('print -dpdf -loose threemonthsonIS26SUT%i',ihc);
%     eval(fileprintpngcmd)
idealresponse
end
