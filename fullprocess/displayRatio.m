%========================== displayratio.m =========================
clear

%=============================================
addpath ZZtoolbox/00gabrielson/
addpath ZZtoolbox/00benoit/
%=======
addpath ZZtoolbox/
%=============================================

% used functions :
%  - smoothpolyLL
%  - theoreticalNRSonIS26.m
%  - 

%=======
for ihc = 1:1
    figure(ihc)
    comload = sprintf('load AAresultswithFB/resultssta26sensor%i',ihc);
    eval(comload);
    
    Dstart = str2double(fileswithdotmat(1).name(13:15));
    Dend   = str2double(fileswithdotmat(length(fileswithdotmat)).name(13:15));
    Dstart = str2double(fileswithdotmat(1).name(13:15));
    Dend   = str2double(fileswithdotmat(length(fileswithdotmat)).name(13:15));
    doubledaynumber = nbmats;
    
    %     permutenbmats = randperm(nbmats);
    %     allRatioPfilters = allRatioPfilters(:,permutenbmats(1:doubledaynumber));
    
    switch ihc
        case 1
            coeffsens=1.037;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 2
            coeffsens=1.065;%.1;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
        case 3
            coeffsens=1.045;
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
    
    allRatioPfilters_ave          = nanmean(allRatioPfilters,2);
    
    
    [allfrqsPfiltersU, inda] = unique(allfrqsPfilters);
    allRatioPfiltersU        = allRatioPfilters(inda,:);
    allmeanMSCcstPfiltersU   = allmeanMSCcstPfilters(inda,:);
    
    allRatioPfiltersUZ        = allRatioPfiltersU(not(allfrqsPfiltersU==0),:);
    allmeanMSCcstPfiltersUZ   = allmeanMSCcstPfiltersU(not(allfrqsPfiltersU==0),:);
    allfrqsPfiltersUZ         = allfrqsPfiltersU(not(allfrqsPfiltersU==0));

    
    if 1
        %======= fit the curve
        segmentsnumber = 2;
        polydegrees    = (0:1:4);
        logtrain.flag  = 1;
        logtrain.N     = 50;
        logfit.flag    = 1;
        logfit.N       = 30;
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
    
%     (xi,yi,P,degs,logtrain,logfit)
    
%     [FreqFitphase_Hz,Rfit] = smoothpolyLL(...
%         allfrqsPfilters(reducefreqrange), ...
%         allRatioPfilters_ave(reducefreqrange),...
%         segmentsnumber,polydegrees,logtrain,logfit);
    [FreqFitphase_Hz,absRfit] = smoothpolyLL(...
        allfrqsPfilters(reducefreqrange), ...
        abs(allRatioPfilters_ave(reducefreqrange)),...
        segmentsnumber,polydegrees,logtrain,logfit);
%     [FreqFitphase_Hz,angleRfit] = smoothpolyLL(...
%         allfrqsPfilters(reducefreqrange), ...
%         angle(allRatioPfilters_ave(reducefreqrange)),...
%         segmentsnumber,polydegrees,logtrain,logfit);

    %=========================================================
%     coeffsens = theoreticalNRSonIS26(Rfit, FreqFitphase_Hz, allRatioPfilters, ...
%         allfrqsPfilters,ihc,ihc+100);

    %     close(ihc+100)

%     figure(ihc+100)
%     subplot(221)
%     title(sprintf('IS26 - sensor #%i, threshold = %4.2f\nday number = %i',...
%         ihc, MSCthreshold, 2*doubledaynumber),'fontname','times','fontsize',12)
%     HorizontalSize = 22;
%     VerticalSize   = 16;
%     set(gcf,'units','centimeters');
%     set(gcf,'paperunits','centimeters');
%     set(gcf,'PaperType','a3');
%     % %             set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
%     set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
%     set(gcf,'color', [1,1,0.92]);
%     set(gcf, 'InvertHardCopy', 'off');
%     
%     fileprintepscmd = sprintf('print -depsc -loose ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
%     fileeps2pdfcmd = sprintf('!epstopdf ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
%     filermcmd = sprintf('!rm ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
    
    
    %==============================================
    figure(ihc)
    clf
%     subplot(211)
    semilogx(allfrqsPfilters,coeffsens*abs(allRatioPfilters),'.', ...
        'color',0.4*ones(3,1))
    hold on
    semilogx(allfrqsPfiltersUZ, coeffsens*nanmean(abs(allRatioPfiltersUZ),2),'.-b')
     semilogx(FreqFitphase_Hz, coeffsens*absRfit,'m','linew',2)
    hold on
    set(gca,'fontname','times','fontsize',16)
    set(gca,'xlim',[0.005 8])
    set(gca,'ylim',[0.9 1.1])
%     semilogx(FreqFitphase_Hz, abs(Rfit),'m','linew',2)
    hold off
    grid on
    xlabel('frequency - Hz')
    ylabel('Estimated Gain')
    drawnow
    
%     subplot(212)
%     semilogx(allfrqsPfiltersU,[allmeanMSCcstPfiltersU])
%     grid on
%     set(gca,'fontname','times','fontsize',10)
%     set(gca,'xlim',[0.001 8])
%     set(gca,'ylim',[0.95 1])
%     xlabel('frequency - Hz')
%     ylabel('MSC')
    
    HorizontalSize = 22;
    VerticalSize   = 16;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a3');
    set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    set(gcf,'color', [1,1,0.92]);
    set(gcf, 'InvertHardCopy', 'off');

    
    printdirectory  = ' ../slidesITW2015/';
    fileprintepscmd = sprintf('print -depsc -loose %s3zones%i.eps',printdirectory,ihc);
    fileeps2pdfcmd  = sprintf('!epstopdf %s3zones%i.eps',printdirectory,ihc);
    filermcmd       = sprintf('!rm %s3monthsonIS26SUT%i.eps',printdirectory,ihc);
    %
%       eval(fileprintepscmd)
%         eval(fileeps2pdfcmd)
    %     eval(filermcmd)
    
    %     fileprintpngcmd = sprintf('print -dpng -loose ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
    %     fileprintpngcmd = sprintf('print -dpdf -loose threemonthsonIS26SUT%i',ihc);
    %     eval(fileprintpngcmd)
    
%     %=====
%     subplot(212)
%     semilogx(allfrqsPfilters,180*angle(allRatioPfilters)/pi,'.', ...
%         'color',0.9*ones(3,1))
%     hold on
% %     semilogx(FreqFitphase_Hz, 180*angleRfit/pi,'m','linew',2)
%     semilogx(allfrqsPfilters, 180*angle(allRatioPfilters_ave)/pi,'r')
% %=======
% %     semilogx(FreqFitphase_Hz, 180*angle(Rfit)/pi,'m','linew',2)
%     hold off
%     set(gca,'fontname','times','fontsize',10)
%     set(gca,'xlim',[0.002 8])
%     set(gca,'ylim',[-10 10])
%     grid on
%     xlabel('frequency - Hz')
%     ylabel('Estimated Phase - degree')
%     drawnow    

    

end
