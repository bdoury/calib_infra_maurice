%========================== displayration.m =========================
clear
% close all
addpath ZZtoolbox/

for ihc = 6:6
    figure(ihc+1)
    comload = sprintf('load AAresults/resultssta26sensor%i',ihc);
    eval(comload);
    Dstart = str2double(filesmat(1).name(13:15));
    Dend   = str2double(filesmat(length(filesmat)).name(13:15));
    nbmats
    doubledaynumber = 10;
    permutenbmats = randperm(nbmats);
    
    allRatioPfilters = allRatioPfilters(:,permutenbmats(1:doubledaynumber));
    
    allRatioPfilters_ave          = nanmean(allRatioPfilters,2);
    absallRatioPfilters_ave       = abs(allRatioPfilters_ave);
    angleallRatioPfilters_ave_deg = angle(allRatioPfilters_ave)*180/pi;

    %======= fit the curve
    segmentsnumber = 2;
    polydegrees    = (0:1:10);
    logtrain.flag  = 1;
    logtrain.N     = 800;
    logfit.flag    = 1;
    logfit.N       = 200;

    [FreqFitabs_Hz,absRfit] = smoothpolyLL(allfrqsPfilters,absallRatioPfilters_ave, ...
        segmentsnumber,polydegrees,logtrain,logfit);

    [FreqFitphase_Hz,phaseRfit_deg] = smoothpolyLL(allfrqsPfilters,angleallRatioPfilters_ave_deg,...
        segmentsnumber,polydegrees,logtrain,logfit);
    
    %=======
    subplot(2,1,1)
    semilogx(allfrqsPfilters,20*log10(abs(allRatioPfilters)),'.', ...
        'color',0.9*ones(3,1))
    hold on
    semilogx(allfrqsPfilters,20*log10(absallRatioPfilters_ave),'.')
    hold on
    semilogx(FreqFitabs_Hz,20*log10(absRfit),'k')
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
    semilogx(FreqFitphase_Hz,phaseRfit_deg,'k')
    hold off
    set(gca,'fontname','times','fontsize',10)
    set(gca,'xlim',[0.002 4])
    set(gca,'ylim',[-5 5])
    grid on
    xlabel('frequency - Hz')
    ylabel('ratio phase - degree')
    drawnow

    subplot(211)
    title(sprintf('IS26 - sensor #%i, threshold = %4.2f\nday number = %i',...
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
    
    
    fileprintepscmd = sprintf('print -depsc -loose ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
    fileeps2pdfcmd = sprintf('!epstopdf ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
    filermcmd = sprintf('!rm ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);

%     eval(fileprintepscmd)
%     eval(fileeps2pdfcmd)
%     eval(filermcmd)

%     fileprintpngcmd = sprintf('print -dpng -loose ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
%     eval(fileprintpngcmd)

end
