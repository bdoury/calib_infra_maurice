clear

addpath ZZtoolbox/

saveflag = 1;
dataresults='AAresultswithFBbis/';
for indexofSTA = [1,3,4,5]
    
            switch indexofSTA
            case 1
                %                         coeffsens=1.048;
                coeffsens=1.015;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 2
                coeffsens=1.1;%.1;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 3
                coeffsens=1.04;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 4
                coeffsens=1.05;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 5
                coeffsens=1.06;
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
    
    
    
    
    filesindir = dir(sprintf('%ss%i/*.mat',dataresults,indexofSTA));
    nbofcouplesdays = length(filesindir);
    allsRatio = zeros(10000,nbofcouplesdays);
    for inday=1:nbofcouplesdays
        filename = filesindir(inday).name;
        comload = sprintf('load %ss%i/%s',dataresults,indexofSTA,filename);
        eval(comload);
        if inday==1
            frqs = allfrqsPfilters;
            allsRatio = allsRatio(1:length(frqs),:);
        end
        allsRatio(:,inday) = allRatioPfilters(:,inday);
    end
    allsRatio             = allsRatio(1:length(frqs),:);
    [frqsU, indunique]    = unique(frqs);
    allRatioU             = allsRatio(indunique,:);
    frqsUZ               = frqsU(not(frqsU==0));
    allRatioUZ            = allRatioU(not(frqsU==0),:);
    meanAllratioUZ        = nanmean(allRatioUZ,2);
    absmeanAllratioUZ     = abs(meanAllratioUZ);
    phasemeanAllratioUZ  = angle(meanAllratioUZ);
    
    segmentsnumber = 2;
    polydegrees    = (0:1:8);
    logtrain.flag  = 1;
    logtrain.N     = 70;
    logfit.flag    = 1;
    logfit.N       = 50;
    
    [FreqFitabs_Hz,absRfit] = smoothpolyLL(frqsU, ...
        abs(nanmean(allRatioU,2)),...
        segmentsnumber,polydegrees,logtrain,logfit);
    
    
    % figure(1)
    %
    % Npermut=10;
    % P = randperm(nbofcouplesdays,Npermut);
    % meanP = abs(nanmean(allRatioU(:,P),2));
    % Q = setdiff((1:nbofcouplesdays),P);
    %     semilogx(frqsU,20*log10(meanP),'r')
    % hold on
    % semilogx(frqsU,20*log10(meanP*1.05),'r')
    % semilogx(frqsU,20*log10(meanP*0.95),'r')
    %
    % for inday=1:nbofcouplesdays-Npermut
    %     semilogx(frqsU,20*log10(abs(allRatioU(:,Q(inday)))),'.','color',0.5*[1 1 1])
    % end
    % hold off
    % grid on
    coeffsens_dB = +20*log10(coeffsens);
    
    figure(indexofSTA)
    PP = randperm(nbofcouplesdays);
    allRatioUZP = allRatioUZ(:,PP);
    
    Nfollowing = 8;
    % semilogx(frqsUZ,20*log10(absmeanAllratioUZ),'r')
    % hold on
    semilogx(frqsUZ,20*log10(absmeanAllratioUZ*1.05)+coeffsens_dB,'r','linew',2)
    hold on
    semilogx(frqsUZ,20*log10(absmeanAllratioUZ*0.95)+coeffsens_dB,'r','linew',2)
    
    % for ii=1:round(nbofcouplesdays/Nfollowing)
    %     id1=(ii-1)*Nfollowing+1;
    %     id2=min([id1+Nfollowing,nbofcouplesdays]) ;
    %     semilogx(frqsUZ,20*log10(abs(nanmean(allRatioUZP(:,id1:id2),2))),'.','color',0.5*[1 1 1])
    %     hold on
    % end
    
    nbdraw = 100;
    for ii=1:nbdraw
        PP_ii = randperm(nbofcouplesdays);
        semilogx(frqsUZ,20*log10(abs(nanmean(allRatioUZP(:,PP_ii(1:Nfollowing)),2)))+coeffsens_dB,'.','color',0.3*[1 1 1])
        hold on
    end
    
    
    
    semilogx(FreqFitabs_Hz, 20*log10(absRfit)+coeffsens_dB,'b')
    hold off
    grid on
    
    set(gca,'fontname','times','fontsize',14)
    set(gca,'xlim',[0.005 10])
    set(gca,'ylim',[-2 2])
    xlabel('Frequency [Hz]')
    ylabel('Gain ratio [dB]')
    title(sprintf('IS26 -  sensor H%i\ndashed line: +/-5%s on the averaged amplitude\n %i points computed on %i random selected days', indexofSTA,'%',nbdraw,2*Nfollowing),...
        'fontname','times','fontsize',14)
    
    
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
    fileprintepscmd = sprintf('print -depsc -loose %srandondrawoncoupleofdays%i.eps',printdirectory,indexofSTA);
    fileeps2pdfcmd  = sprintf('!epstopdf %srandondrawoncoupleofdays%i.eps',printdirectory,indexofSTA);
    filermcmd       = sprintf('!rm %srandondrawoncoupleofdays%i.eps',printdirectory,indexofSTA);
    
    if saveflag
        eval(fileprintepscmd)
        eval(fileeps2pdfcmd)
        eval(filermcmd)
    end
end
