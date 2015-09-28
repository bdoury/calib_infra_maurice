clear
close all
for ihc=1:5
    comload = sprintf('load AAresultsBW0812Hz/resultssta26sensor%i.mat',ihc);
    eval(comload)
    
    Dstart = str2double(fileswithdotmat(1).name(13:15));
    Dend   = str2double(fileswithdotmat(length(fileswithdotmat)).name(13:15));
    if 1
        doubledaynumber = (Dend-Dstart+3)/2;
    else
        doubledaynumber = 10;
        permutenbmats = randperm(nbmats);
        allRatioPfilters = allRatioPfilters(:,permutenbmats(1:doubledaynumber));
    end
    
    indselect = find(allfrqsPfilters(:,1)>=1,1,'first');
    allfrqsPfilters(indselect)
    figure(ihc)
    subplot(211);
    plot(abs(allRatioPfilters(indselect,:)),'or','markersize',4)
    set(gca,'ylim',[0.7 1.3])
    grid on
    hold on
    plot(ones(2,1)*(1:size(allRatioPfilters,2)),[abs(allRatioPfilters(indselect,:))-allSTDmodRatioPfilters(indselect,:);...)
        abs(allRatioPfilters(indselect,:))+allSTDmodRatioPfilters(indselect,:)],'o-b');
    hold off
    title(sprintf('IS26 -  sensor #%i, threshold = %4.2f\nday number = %i',...
        ihc, MSCthreshold, 2*doubledaynumber),'fontname','times','fontsize',12)
    
    subplot(212);
    plot(allmeanMSCcstPfilters(indselect,:),'ob','markersize',6,...
        'markerfacec','b')
    ylabel('MSC')
    grid on
    
    
    HorizontalSize = 16;
    VerticalSize   = 12;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a3');
    set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    set(gcf,'color', [1,1,0.92]);
    set(gcf, 'InvertHardCopy', 'off');
    
    printdirectory = ' ../slidesITW2015/';
    fileprintepscmd = sprintf('print -depsc -loose %sevolutionon%iatfreq%i.eps',printdirectory,ihc,allfrqsPfilters(indselect));
    
%     eval(fileprintepscmd)
    
end
