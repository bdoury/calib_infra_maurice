%============= temporalevolution.m ============
% ratio at different frequency values,
% The program displays either the sequence of
% the pairs of days or the randomly selection
% of a few number of them. We also apply
% trimmed mean.
%=====================================================
clear
addpath ../ZZtoolbox/
addpath ../ZZtoolbox/00gabrielson
addpath ../ZZtoolbox/00benoit

directorydata = '../AAresultswithFB98/';
printdirectory  = ' ../../figures/';
saveflag = 0;
bootdraw = 1;
nbdraw   = 100;
Ncouples = 15;
trimmingpercent = 30;
numberthreshold = 0;

cp=0;
for  selectedfrequency_Hz = [0.05 0.1 0.5 1]
    cp=cp+1;
    for ihc=2
        comload = sprintf('load %sresultssta26sensor%i.mat',directorydata,ihc);
        numfig = ihc;
        
        eval(comload)
        freqindselect   = find(...
            allfrqsPfilters(:,1)>=selectedfrequency_Hz,1,'first');
        
        if ihc==4
            allpraticalvaluesatselectfreq_Hz = ...
                abs(allRatioSupPfilters(freqindselect,[1:65 67:nbmats]));
        else
            allpraticalvaluesatselectfreq_Hz = ...
                abs(allRatioSupPfilters(freqindselect,:));
        end
        
        nbofcouplesdays = length(allpraticalvaluesatselectfreq_Hz);
        %===================== plots ==============
        figure(numfig)
        if bootdraw
            subplot(2,2,cp)
            doubledaynumber = Ncouples;
            meanallpraticalvaluesat1Hz_ii=zeros(nbdraw,1);
            for idraw=1:nbdraw
                PP_ii = randperm(nbofcouplesdays);
                meanallpraticalvaluesat1Hz_ii(idraw)= ...
                    trimmean(allpraticalvaluesatselectfreq_Hz(:,PP_ii(1:Ncouples)),trimmingpercent);
            end
            mmcenter = nanmean(meanallpraticalvaluesat1Hz_ii);
            allpraticalvaluesatselectfreq_Hz = ...
                meanallpraticalvaluesat1Hz_ii-mmcenter+1;
        else
            doubledaynumber = nbmats;
            subplot(211);
            indnotsufficient = ...
                find(nbofvaluesoverthreshold(freqindselect,:)<=numberthreshold);
            allpraticalvaluesatselectfreq_Hz(indnotsufficient) = NaN;
            mmcenter = nanmean(allpraticalvaluesatselectfreq_Hz);
            allpraticalvaluesatselectfreq_Hz = ...
                allpraticalvaluesatselectfreq_Hz-mmcenter+1;
        end
        plot((1:nbdraw),20*log10(allpraticalvaluesatselectfreq_Hz),...
            'ok','markersize',6,'markerfacec','k')
        hold on
        plot([0 10000], (20*log10(0.95))*[1 1],'--r','linew',2)
        plot([0 10000], (20*log10(1.05))*[1 1],'--r','linew',2)
        hold off
        %     set(gca,'fontname','times','fontsize',14)
        ylabel(sprintf('Gain at %4.2f Hz',selectedfrequency_Hz),...
            'fontname','times','fontsize',12)
        grid on
%         xlabel('days','fontname','times','fontsize',12)
        xlabel('draw number','fontname','times','fontsize',12)
        set(gca,'fontname','times','fontsize',12)
        set(gca,'xlim',[1 nbdraw])
        set(gca,'ylim',[-1 1])      
%         title(sprintf('IS26 -  sensor #%i, day number = %i, threshold = %4.2f',...
%             ihc, 2*doubledaynumber, MSCthreshold),'fontname','times','fontsize',12)
        
        if cp==1
            title(sprintf('station %i, day number = %i\ntrimming at %i%s',...
                ihc, 2*doubledaynumber,trimmingpercent,'%'),'fontname','times','fontsize',12)
        end
        hold off
        if not(bootdraw)
            subplot(212);
            %             plot(allmeanMSCcstPfilters(indselect,:),'ob','markersize',6,...
            %                 'markerfacec','b')
            semilogy((1:nbmats), nbofvaluesoverthreshold(freqindselect,:),'ob','markersize',6,...
                'markerfacec','b')
            %             set(gca,'ylim',[0.98 1])
            ylabel('number fo values','fontname','times','fontsize',14)
            grid on
            set(gca,'fontname','times','fontsize',12)
            set(gca,'xlim',[1 2*nbofcouplesdays])
            xlabel('days','fontname','times','fontsize',12)
        end
        hold off
    end
end
HorizontalSize = 16;
VerticalSize   = 14;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');

fileprint = sprintf('%sevolutionon%iatdifffreq.eps',printdirectory,ihc);
% fileprint = sprintf('%sevolutionon%iatfreq1Hz.eps',printdirectory,ihc);

fileprintepscmd = sprintf('print -depsc -loose %s',fileprint);
fileeps2pdfcmd  = sprintf('!epstopdf %s',fileprint);
filermcmd       = sprintf('!rm %s',fileprint);
%
if saveflag
    eval(fileprintepscmd)
    eval(fileeps2pdfcmd)
    eval(filermcmd)
end

