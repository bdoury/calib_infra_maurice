%========================== displayratio.m =========================
clear

%=============================================
addpath ZZtoolbox/00gabrielson/
<<<<<<< HEAD
addpath ZZtoolbox/00benoit/
=======
>>>>>>> e518a59897f74aa44d579d6b4e37dafce4bba1aa
addpath ZZtoolbox/
%=============================================

% used functions :
%  - smoothpolyLL
%  - theoreticalNRSonIS26.m
%  - 

<<<<<<< HEAD
for ihc = 8:8
=======
for ihc = 6:8
>>>>>>> e518a59897f74aa44d579d6b4e37dafce4bba1aa
    figure(ihc)
    comload = sprintf('load AAresults/resultssta26sensor%i',ihc);
    eval(comload);
    
    
<<<<<<< HEAD
    Dstart = str2double(filesmat(1).name(13:15));
    Dend   = str2double(filesmat(length(filesmat)).name(13:15));
=======
    Dstart = str2double(fileswithdotmat(1).name(13:15));
    Dend   = str2double(fileswithdotmat(length(fileswithdotmat)).name(13:15));
>>>>>>> e518a59897f74aa44d579d6b4e37dafce4bba1aa
    doubledaynumber = nbmats;
    
    %     permutenbmats = randperm(nbmats);
    %     allRatioPfilters = allRatioPfilters(:,permutenbmats(1:doubledaynumber));
    
    allRatioPfilters_ave          = nanmean(allRatioPfilters,2);
    
    if 1
        %======= fit the curve
        segmentsnumber = 2;
<<<<<<< HEAD
        polydegrees    = (0:1:7);
        logtrain.flag  = 1;
        logtrain.N     = 500;
        logfit.flag    = 1;
        logfit.N       = 300;
=======
        polydegrees    = (0:1:9);
        logtrain.flag  = 1;
        logtrain.N     = 500;
        logfit.flag    = 1;
        logfit.N       = 200;
>>>>>>> e518a59897f74aa44d579d6b4e37dafce4bba1aa
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
    
    [FreqFitphase_Hz,Rfit] = smoothpolyLL(...
        allfrqsPfilters(reducefreqrange), ...
        allRatioPfilters_ave(reducefreqrange),...
        segmentsnumber,polydegrees,logtrain,logfit);
<<<<<<< HEAD
    [FreqFitphase_Hz,absRfit] = smoothpolyLL(...
        allfrqsPfilters(reducefreqrange), ...
        abs(allRatioPfilters_ave(reducefreqrange)),...
        segmentsnumber,polydegrees,logtrain,logfit);
    [FreqFitphase_Hz,angleRfit] = smoothpolyLL(...
        allfrqsPfilters(reducefreqrange), ...
        angle(allRatioPfilters_ave(reducefreqrange)),...
        segmentsnumber,polydegrees,logtrain,logfit);
=======
>>>>>>> e518a59897f74aa44d579d6b4e37dafce4bba1aa

    %=========================================================
    coeffsens = theoreticalNRSonIS26(Rfit, FreqFitphase_Hz, allRatioPfilters, ...
        allfrqsPfilters,ihc,ihc+100);
<<<<<<< HEAD
%     close(ihc+100)
=======
    
>>>>>>> e518a59897f74aa44d579d6b4e37dafce4bba1aa
    figure(ihc+100)
    subplot(221)
    title(sprintf('IS26 - sensor #%i, threshold = %4.2f\nday number = %i',...
        ihc, MSCthreshold, 2*doubledaynumber),'fontname','times','fontsize',12)
    HorizontalSize = 22;
    VerticalSize   = 16;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a3');
    % %             set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    set(gcf,'color', [1,1,0.92]);
    set(gcf, 'InvertHardCopy', 'off');
    
    fileprintepscmd = sprintf('print -depsc -loose ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
    fileeps2pdfcmd = sprintf('!epstopdf ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
    filermcmd = sprintf('!rm ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
    
    
    %==============================================
    figure(ihc)
    subplot(211)
    semilogx(allfrqsPfilters,abs(allRatioPfilters),'.', ...
        'color',0.9*ones(3,1))
    hold on
<<<<<<< HEAD
    semilogx(allfrqsPfilters, abs(allRatioPfilters_ave),'.r')
    semilogx(FreqFitphase_Hz, absRfit,'m','linew',2)
    hold off
    set(gca,'fontname','times','fontsize',10)
    set(gca,'xlim',[0.002 8])
    set(gca,'ylim',[0.8 1.1])
=======
    semilogx(FreqFitphase_Hz, abs(Rfit),'m','linew',2)
    hold off
    set(gca,'fontname','times','fontsize',10)
    set(gca,'xlim',[0.002 8])
    set(gca,'ylim',[0.9 1.1])
>>>>>>> e518a59897f74aa44d579d6b4e37dafce4bba1aa
    grid on
    xlabel('frequency - Hz')
    ylabel('Estimated Gain')
    drawnow
    %=====
    subplot(212)
    semilogx(allfrqsPfilters,180*angle(allRatioPfilters)/pi,'.', ...
        'color',0.9*ones(3,1))
    hold on
<<<<<<< HEAD
    semilogx(FreqFitphase_Hz, 180*angleRfit/pi,'m','linew',2)
    semilogx(allfrqsPfilters, 180*angle(allRatioPfilters_ave)/pi,'r')
=======
    semilogx(FreqFitphase_Hz, 180*angle(Rfit)/pi,'m','linew',2)
>>>>>>> e518a59897f74aa44d579d6b4e37dafce4bba1aa
    hold off
    set(gca,'fontname','times','fontsize',10)
    set(gca,'xlim',[0.002 8])
    set(gca,'ylim',[-10 10])
    grid on
    xlabel('frequency - Hz')
    ylabel('Estimated Phase - degree')
    drawnow

end
