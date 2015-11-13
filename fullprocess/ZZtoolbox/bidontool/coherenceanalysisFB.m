clear
allcols = ['b';'k';'r';'c';'m';'g';'k';'r';'c';'m';'g';'k';'r';'c';'m';'g'];
for itour=1:1
    col=allcols(itour,:);
    highB = 0;
    lowB = 1;
    addpath  ../../../ZZtoolbox/
    %=====================
    % under test = 1, reference = 2
    %=====================
    %
    if 1
        %===================== read data =========================
        directorydata = '../../../../DATA_IS/';
        % directorydata = './';
        filesmat = dir([directorydata 'I26*.mat']);
        nbmats   = length(filesmat);
        id1=1;
        for ifile=1:1%nbmats
            file_i=filesmat(ifile).name;
            dotpos = strfind(file_i,'.');
            underscorepos = strfind(file_i,'_');
            fileonly = file_i(setdiff(1:dotpos-1,underscorepos));
            commandload    = sprintf('load %s%s',directorydata,file_i);
            eval(commandload)
            %     Ts_sec = 1/Fs_Hz;
            if exist('data_pa','var')
                signals        = data_pa;
                Ts_sec         = (time_pa(2)-time_pa(1))*24*3600;
            elseif exist('data','var'),
                signals = data;
                Ts_sec         = (time(2)-time(1))*24*3600;
            elseif exist('records','var'),
                Lrecords=length(records)/2;
                Ts_sec  = 1/samprate(1);
                for idchan = 1:Lrecords
                    % if L is the number of records,
                    % the L/2 first records are H (under test)
                    % and the others (add L/2) are C (reference)
                    LL = min([length(records{idchan}.data), ...
                        length(records{idchan+Lrecords}.data)]);
                    id2=id1+LL-1;
                    signals(id1:id2,2) = records{idchan}.data(1:LL);
                    signals(id1:id2,1) = records{idchan+Lrecords}.data(1:LL);
                    id1=id2+1;
                end
            end
        end
    else
        clear signals
        ht = [1 1]/2;
        T = 1000000;
        e = randn(T,1);
        sigmaUT = 0.05;
        sigmaREF = sqrt(60)* sigmaUT;
        eUT = e+sigmaUT*randn(T,1);
        eREF = e+sigmaREF*randn(T,1);
        signals(:,1)=filter(ht,1,eUT);
        signals(:,2)=filter(ht,1,eREF);
        Ts_sec = 1/20;
    end
    Fs_Hz     = 1/Ts_sec;
    %============================================
    
    %============================================
    beginband = 0;
    %====
%     if exist('ht','var')
%          HUT_true = fft(ht,Lfft_ii);
%         SUU_true = (abs(HUT_true) .^2)*(1+sigmaUT^2);
%         SRR_true = (1+sigmaREF^2);
%         SUR_true = HUT_true;
%         MSC_true = abs(SUR_true) .^2 ./ (SRR_true .* SUU_true);
%     end
    %===================== initial values ======================
    signals = double(signals);
    T = size(signals,1);
    Nlogband = 3;
    signalsH = signals;
    signalsL = signals;
    Tstationarity_sec   = 1000;
    TaverageFFT_sec_ii  = Tstationarity_sec / 2^(Nlogband);
    overlapFFT     = 0.5;
    overlapAVE     = 0.0;
    Fs_Hz_ii = Fs_Hz;
    %===================== filter bank =========================
    [LO_Decomp,HI_Decomp,LO_Recons,HI_Recons] = wfilters('db32');
    for iband = 1:Nlogband
        signalsH  = filter(HI_Decomp,1,signalsL);
        % decimation
        signalsH  = signalsH(1:2:end,:);
        [T_ii,nbsignals]  = size(signalsH);
        signalsH  = signalsH-ones(T_ii,1)*mean(signalsH);
        
        signalsL  = filter(LO_Decomp,1,signalsL);
        signalsL  = signalsL(1:2:end,:);
        signalsL  = signalsL-ones(T_ii,1)*mean(signalsL);
        
        % we divide by 2
        Fs_Hz_ii    = Fs_Hz_ii/2;
        TaverageFFT_sec_ii = TaverageFFT_sec_ii*2;TaverageFFT_sec_ii*Fs_Hz_ii
        Tfft_sec_ii  = TaverageFFT_sec_ii/20;
        
        Lfft_ii      = 2*round(Tfft_sec_ii*Fs_Hz_ii);

        %
        %===================== spectral analysis ===================

        [allSDs, time_sec_ii, frqsFFT_Hz_ii] = ...
            estimSDs(signalsH(:,1),signalsH(:,2),...
            ones(Lfft_ii,1),overlapFFT, ...
            TaverageFFT_sec_ii, overlapAVE, Fs_Hz_ii, 'hann');
        nbblocksAVE = length(allSDs);
        frqsFFT_Hz_ii_aff=frqsFFT_Hz_ii+Fs_Hz_ii/2;
        
        %%
        %================================================================
        coherenceThreshold = 0.99;
        %================================================================
        indextabMSCsupthreshold = false(Lfft_ii,nbblocksAVE);
        indextabMSCsupthreshold([allSDs.MSC]>coherenceThreshold) = true;
        
        tabSUU = [allSDs.UU];
        tabSUUwithMSCsupeta = nan(Lfft_ii,nbblocksAVE);
        tabSUUwithMSCsupeta(indextabMSCsupthreshold) = ...
            tabSUU(indextabMSCsupthreshold);
        SUUmedian = nanmean(tabSUUwithMSCsupeta,2);
        
        tabSRR = [allSDs.RR];
        tabSRRwithMSCsupeta = nan(Lfft_ii,nbblocksAVE);
        tabSRRwithMSCsupeta(indextabMSCsupthreshold) = ...
            tabSRR(indextabMSCsupthreshold);
        SRRmedian = nanmean(tabSRRwithMSCsupeta,2);
        
        tabSUR = [allSDs.UR];
        tabSURwithMSCsupeta = nan(Lfft_ii,nbblocksAVE);
        tabSURwithMSCsupeta(indextabMSCsupthreshold) = ...
            tabSUR(indextabMSCsupthreshold);
        SURmedian = nanmean(tabSURwithMSCsupeta,2);
        
        tabMSC = [allSDs.MSC];
        tabMSCwithMSCsupeta = nan(Lfft_ii,nbblocksAVE);
        tabMSCwithMSCsupeta(indextabMSCsupthreshold) = ...
            tabMSC(indextabMSCsupthreshold);
        
        allT.TUUonUR                         = linspace(0.6,3,100);
        allT.TURonRR                         = linspace(0.6,3,100);
        allT.MSC                             = linspace(0.5,1,100);
        %================= on HUUUR
        tabHUUUR = [allSDs.Rsup];
        tabmodHUUUR = (abs(tabHUUUR));
        hatHUUUR = nanmedian(tabmodHUUUR,2);
        
        tabmodHUUURwithMSCsupeta = nan(Lfft_ii,nbblocksAVE);
        tabmodHUUURwithMSCsupeta(indextabMSCsupthreshold) = ...
            (abs(tabHUUUR(indextabMSCsupthreshold)));
        hatHUUURwithMSCsupeta = nanmedian(tabmodHUUURwithMSCsupeta,2);
        
        tabHURRR = [allSDs.Rinf];
        tabmodHURRRwithMSCsupeta = nan(Lfft_ii,nbblocksAVE);
        tabmodHURRRwithMSCsupeta(indextabMSCsupthreshold) = ...
            (abs(tabHURRR(indextabMSCsupthreshold)));
        
        hatHURRRwithMSCsupeta = nanmedian(tabmodHURRRwithMSCsupeta,2);
        
        %================= on HURRR
        tabHURRR = [allSDs.Rinf];
        tabmodHURRRwithMSCsupeta = nan(Lfft_ii,nbblocksAVE);
        tabmodHURRRwithMSCsupeta(indextabMSCsupthreshold) = ...
            (abs(tabHURRR(indextabMSCsupthreshold)));
        
        %     %%itour
        %     figure(1)
        %     clf
        %     semilogx(frqsFFT_Hz,20*log10(hatHUUUR),'k')
        %     hold on
        %     semilogx(frqsFFT_Hz,20*log10(hatHUUURwithMSCsupeta),'r')
        %     set(gca,'xlim',[frqsFFT_Hz(2) 6])
        %     set(gca,'ylim',[-6 2])
        %     hold off
        %     grid on
        %
        %
        %     %%
        %     Gsup as Ginf
        %     figure(7)
        %     subplot(211)
        %     semilogx(frqsFFT_Hz,20*log10(hatHUUUR),'k')
        %     hold on
        %     semilogx(frqsFFT_Hz,20*log10(hatHURRR),'r')
        %     set(gca,'xlim',[frqsFFT_Hz(2) 6])
        %     set(gca,'ylim',[-6 2])
        %     hold off
        %     grid on
        %     hold on
        %     hdcenterfreqFFTH = semilogx(...
        %         sum(frqsFFT_Hz(2:3))*ones(2,1)/2,[-1 1]);
        %     hold off
        %
        %     for ib=2:fix(fix(Lfft_ii/2))
        %         subplot(212)
        %         %         hist(tabmodHURRRwithMSCsupeta(ib,:),50)
        %         plot(tabmodHURRRwithMSCsupeta(ib,:),tabmodHUUURwithMSCsupeta(ib,:),'.')
        %         hold on
        %         plot(median(tabmodHURRRwithMSCsupeta(ib,:)), ...
        %             median(tabmodHUUURwithMSCsupeta(ib,:)),'or',...
        %             'markerfac','y','markers',10)
        %         plot(1, 1,'ok','markerfac','k','markers',10)
        %         %         plot(mean(tabmodHURRRwithMSCsupeta(ib,:)), ...
        %         %             mean(tabmodHUUURwithMSCsupeta(ib,:)),'ok',...
        %         %             'markerfac','k','markers',6)
        %         set(gca,'xlim',[0.6 1.1],'ylim',[0.8 1.2])
        %         grid on
        %         hold off
        %         drawnow
        %         subplot(211)
        %         if ib<Lfft-1
        %             delete(hdcenterfreqFFTH)
        %             hold on
        %             hdcenterfreqFFTH = semilogx(...
        %                 sum(frqsFFT_Hz(ib:ib+1))*ones(2,1)/2,[-1 1]);
        %
        %             hold off
        %         end
        %         pause
        %         subplot(212)
        %     end
        %%
        % %======= boxplot of Gsup and hist of MSC
        %     figure(8)
        %     subplot(211)
        %     semilogx(frqsFFT_Hz,20*log10(hatHUUUR),'k')
        %     set(gca,'xlim',[frqsFFT_Hz(2) 6])
        %     set(gca,'ylim',[-2 2])
        %     hold off
        %     grid on
        %     hold on
        %     hdcenterfreqFFTL = semilogx(...
        %         mean(frqsFFT_Hz(1:2))*ones(2,1),[-2 2]);
        %     hdcenterfreqFFTH = semilogx(...
        %         mean(frqsFFT_Hz(2:3))*ones(2,1),[-2 2]);
        %     hold off
        %
        %     for ib=2:fix(fix(Lfft/2))
        %         subplot(223)
        %         hist(tabMSCwithMSCsupeta(ib,:))
        %         drawnow
        %         subplot(224)
        %         boxplot(tabmodHUUURwithMSCsupeta(ib,:))
        %         drawnow
        %         set(gca,'ylim',[0 100])
        %         subplot(211)
        %         if ib<Lfft-1
        %             delete(hdcenterfreqFFTH)
        %             delete(hdcenterfreqFFTL)
        %             hold on
        %             hdcenterfreqFFTL = semilogx(...
        %                 mean(frqsFFT_Hz(ib-1:ib))*ones(2,1),[-2 2]);
        %             hdcenterfreqFFTH = semilogx(...
        %                 mean(frqsFFT_Hz(ib:ib+1))*ones(2,1),[-2 2]);
        %             hold off
        %         end
        % %         pause
        %         subplot(212)
        %     end
        %%
        % %=======  hist of MSC
        %     figure(9)
        %     clf
        %     subplot(211)
        %     semilogx(frqsFFT_Hz,20*log10(hatHUUUR),'k')
        %     hold on
        %     semilogx(frqsFFT_Hz,20*log10(hatHUUURwithMSCsupeta),'r')
        %     set(gca,'xlim',[frqsFFT_Hz(2) 6])
        %     set(gca,'ylim',[-6 2])
        %     hold off
        %     grid on
        %     hold on
        %     hdcenterfreqFFTL = semilogx(...
        %         mean(frqsFFT_Hz(1:2))*ones(2,1),100*[-1 1]);
        %     hdcenterfreqFFTH = semilogx(...
        %         mean(frqsFFT_Hz(2:3))*ones(2,1),100*[-1 1]);
        %     hold off
        %     rangeMSC=linspace(0.1,0.99999,100);
        %     for ib=2:Lfft/4
        %         subplot(212)
        %         [hbinMSC, binMSC] = hist((tabMSC(ib,:)),20);
        %         hatpdfMSC = hbinMSC/nbblocksAVE/(binMSC(2)-binMSC(1));
        %         bar(binMSC,hatpdfMSC)
        %          pdfMSCth = pdfMSC(rangeMSC',...
        %              nanmedian(tabMSC(ib,:)),...
        %              NaverageFFTs);
        % %         pdfMSCth = pdfMSC(binMSC',nanmean(MSC_true),NaverageFFTs/2);
        %         hold on
        %         plot(rangeMSC,pdfMSCth,'.-r')
        %         hold off
        %         drawnow
        %
        %         subplot(211)
        %         if ib<Lfft-1
        %             delete(hdcenterfreqFFTH)
        %             delete(hdcenterfreqFFTL)
        %             hold on
        %             hdcenterfreqFFTL = semilogx(...
        %                 mean(frqsFFT_Hz(ib-1:ib))*ones(2,1),100*[-1 1]);
        %             hdcenterfreqFFTH = semilogx(...
        %                 mean(frqsFFT_Hz(ib:ib+1))*ones(2,1),100*[-1 1]);
        %             hold off
        %         end
        %         drawnow
        % %         pause
        %         subplot(212)
        %     end
        %%
        %     figure (7)
        %     hold on
        %     SUUtmp = allSDs.UU;
        %     SRRtmp = allSDs.RR;
        %     plot(frqsFFT_Hz(2:fix(fix(Lfft/2)))+beginband,mean(10*log10(SUUtmp(2:fix(fix(Lfft/2)))),2))
        %     hold on
        %     plot(frqsFFT_Hz(2:fix(fix(Lfft/2)))+beginband,mean(10*log10(SRRtmp(2:fix(fix(Lfft/2)))),2))
        %     hold off
        %     set(gca,'xlim',[0 10])
        %     set(gca,'ylim',[10 80])
        %     grid on
        %%
        %    figure (8)
        %     clf
        %     % plot(tabMSCwithMSCsupeta(10,:),tabmodHUUURwithMSCsupeta(10,:),'.')
        %     plot(tabMSCwithMSCsupeta(2:1000,:)',20*log10(tabmodHUUURwithMSCsupeta(2:1000,:)'),'.')
        %     %     boxplot(tabmodHURRRwithMSCsupeta(2:10:1000,:)')
        %     boxplot(abs(tabMSCwithMSCsupeta(2:100:fix(fix(fix(Lfft/2))),:)'))
        %
        %     %     semilogx(frqsFFT_Hz,nanmean(tabMSC,2))
        %     %     hold on
        %     %     semilogx(frqsFFT_Hz,MSC_true,'m')
        %     %     hold off
        %
        %     set(gca,'ylim',[0 1.2])
        %     set(gca,'YScale','lin')
        %     set(gca,'XScale','log')
        %     grid on
        %     %%
        %     figure (10)
        % %     hold on
        %     semilogx(frqsFFT_Hz,20*log10(abs(hatHUUUR)),'.-k')
        % %         hold on
        % %         semilogx(frqsFFT_Hz,20*log10(abs(hatHURRR)),'.-r')
        % %         hold off
        %     set(gca,'xlim',[0 Fs_Hz/2])
        %     set(gca,'ylim',[-6 3])
        %     set(gca,'YScale','lin')
        % %     set(gca,'XScale','log')
        %     grid on
        % %
        %     subplot(313)
        %     hb=boxplot(abs(tabmodHUUURwithMSCsupeta(2:1:1000,:))');
        %     set(gca,'ylim',[0.9 1.1])
        %     set(gca,'YScale','lin')
        %     set(gca,'XScale','log')
        %     grid on
        %
        %%
        % %     %%% ================== check the theoretical probability density
        % %         toto=zeros(fix(fix(Lfft/2)),1);
        %         for indfrq = 2:10:fix(fix(Lfft/2)),indfrq
        %             Rib = [SUUmedian(indfrq) SURmedian(indfrq); conj(SURmedian(indfrq)) SRRmedian(indfrq) ];
        %             if abs(det(Rib))>1
        %                 Ribdiagout = 0.99*sqrt(SUUmedian(indfrq)*SRRmedian(indfrq));
        %                 Rib = [SUUmedian(indfrq) Ribdiagout; Ribdiagout SRRmedian(indfrq) ];
        %             end
        %             [statUUonUR, statURonRR, statMSC] = statsRatiosHbis(allT, Rib, (NaverageFFTs)/2, 0.7);
        %             if not(isempty(statUUonUR.median))
        %             toto(indfrq)=statUUonUR.median;
        %             end
        %         end
        %             expHUUUR= tabmodHUUURwithMSCsupeta(indfrq,:);
        %             if 1
        %                 [hHUUUR, binUUUR] = hist(expHUUUR,20);
        %                 hatpdfHUUUR = hHUUUR/length(expHUUUR)/(binUUUR(2)-binUUUR(1));
        %                 bar(binUUUR, hatpdfHUUUR)
        %                 hold on
        %                 plot(allT.TUUonUR,statUUonUR.pdf,'r')
        %                 hold off
        %                 set(gca,'ylim',[0 6])
        %                 set(gca,'xlim',[0.5 1.8])
        %                 drawnow
        %                 frqsFFT_Hz(indfrq)
        %                 pause
        %             else
        %                 expMSC = tabMSC(indfrq,:);
        %                 [hMSC, binMSC] = hist(expMSC,20);
        %                 hatpdfMSC = hMSC/length(expMSC)/(binMSC(2)-binMSC(1));
        %                 bar(binMSC, hatpdfMSC)
        %                 hold on
        %                 plot(allT.MSC,statMSC.pdf,'r')
        %                 hold off
        %                 set(gca,'ylim',[0 12])
        %                 set(gca,'xlim',[0.5 1])
        %                 drawnow
        %                 frqsFFT_Hz(indfrq)
        %                 pause
        %             end
        %
        %%
        
        % %====================================================================
        % indexrangeFit=2:800;
        % indfreqCut_Hz = 0.4;
        % [HmedianSmo, frqsHmedianSmo] = smoothbypoly((Hmedian(indexrangeFit)),frqsFFT_Hz(indexrangeFit),...
        %     indfreqCut_Hz,[5 5]);
        % %==========================
        % indexrangeFit=2:800;
        % indfreqCut_Hz = 0.4;
        % [phaseSmo, frqsphaseSmo] = smoothbypoly((phasemedian(indexrangeFit)),frqsFFT_Hz(indexrangeFit),...
        %     indfreqCut_Hz,[5 5]);
        %%
        %=============================== pcolor of MSC, etc
        figure(16)
%         subplot(311)
        if not(iband==1)
            hold on
        else
            hold off 
        end
        pcolor(frqsFFT_Hz_ii_aff(1:fix(Lfft_ii/2)),[time_sec_ii.SD]/3600,tabMSC(1:fix(Lfft_ii/2),:)');
        ylabel('time - hour','fontsi',13,'fontn','times')
        shading flat
        set(gca,'XScale','lin')
        set(gca,'xlim',[0.003 10])
        set(gca,'xdir','normal')
        grid on
%         title(sprintf('%s - all MSCs',fileonly(1:5)),'fontsi',13)
        set(gca,'fontn','times','fonts',10)
        colorbar
        %         set(gca,'position',hbox211);
        
        %         subplot(212)
        %         boxplot(tabMSC,'plotstyle','compact')
        %             set(gca,'XScale','log')
        %         grid on
        %         set(gca,'xtick',[],'xtickLabel',[])
        
        %         subplot(412)
        %         pcolor(frqsFFT_Hz,[time_sec.SD]/3600,double(indextabMSCsupthreshold)');
        %         caxis([0 2])
        %         shading flat
        %         ylabel('time - hour','fontsi',16,'fontn','times')
        %         xlabel('frequency - Hz','fontsi',16,'fontn','times')
        %         set(gca,'XScaleuntitled.m','log')
        %         set(gca,'xlim',[0 1.05*Fs_Hz/2])
        %         set(gca,'xdir','norm')
        %         grid on
        %         title(sprintf('MSC>%4.2f',coherenceThreshold),'fontsi',16,'fontn','times')
        %         set(gca,'fontn','times','fonts',10,'fontn','times')
        %
        
%         % subplot(413)
%         if not(itour==1)
%             hold on
%         else 
%             hold off
%         end
%         subplot(312)
%         hold on
%         %         hbox212=[ 0.13    0.45    0.67    0.2];
%         semilogx(frqsFFT_Hz_ii(1:fix(Lfft_ii/2)),20*log10(hatHUUURwithMSCsupeta(1:fix(Lfft_ii/2),:)),'r');
%         % hold on
%         % semilogx(frqsHmedianSmo,(20*log10(HmedianSmo)),'.-r')
%         set(gca,'ylim',[-5 1])
%         set(gca,'xlim',[0.003 8])
%         grid on
%         hold off
%         title(sprintf('Gains with median over the time slots with MSC > %4.2f',...
%             coherenceThreshold),...
%             'fontsi',13,'fontn','times')
%         set(gca,'fontn','times','fonts',10,'fontn','times')
%         ylabel('dB','fontsi',13,'fontn','times')
%         xlabel('frequency - Hz','fontsi',13,'fontn','times')
%         %         set(gca,'Position',hbox212)
        %%
        
        %     subplot(414)
        %     semilogx(frqsFFT_Hz,phasemedian*180/pi,'.-')
        %     hold on
        %     semilogx(frqsHmedianSmo,phaseSmo*180/pi,'.-r')
        %     set(gca,'ylim',[-60 60])
        %     set(gca,'xlim',[0 1.05*Fs_Hz/2])
        %     grid on
        %     hold off
        %     title(sprintf('Phase medians over the time slots with MSC > %4.2f',coherenceThreshold),'fontsi',14,'fontn','times')
        %     set(gca,'fontn','times','fonts',10,'fontn','times')
        %     ylabel('degree','fontsi',16,'fontn','times')
        %     xlabel('frequency - Hz','fontsi',16,'fontn','times')
        
        HorizontalSize = 18;
        VerticalSize   = 25;
        set(gcf,'units','centimeters');
        set(gcf,'paperunits','centimeters');
        set(gcf,'PaperType','a4');
%         set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
        set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
        
        set(gcf,'color', [1,1,0.92]);
        set(gcf, 'InvertHardCopy', 'off');
        
        % figure(8)
        % commandprint = sprintf('print -dpdf %s',fileonly);
        % eval(commandprint)
        
        %%
        % figure(2)
        % clf
        %
        % subplot(211)
        % semilogx(frqsFFT_Hz,20*log10(tabmodHUUURwithMSCsupeta))
        % set(gca,'ylim',[-15 5])
        % set(gca,'xlim',[0 6])
        % grid on
        % ylabel('dB')
        % set(gca,'fontn','times','fonts',10)
        %
        % subplot(212)
        % semilogx(frqsFFT_Hz,tabphaseHUUURwithMSCsupeta)
        % set(gca,'xlim',[0 6])
        % set(gca,'ylim',[-15 15])
        % grid on
        % ylabel('degree')
        % set(gca,'fontn','times','fonts',10)
        %
        %
        % HorizontalSize = 18;
        % VerticalSize   = 18;
        % set(gcf,'units','centimeters');
        % set(gcf,'paperunits','centimeters');
        % set(gcf,'PaperType','a4');
        % set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
        % set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
        %
        % set(gcf,'color', [1,1,0.92]);
        % set(gcf, 'InvertHardCopy', 'off');
        %
        %
        % % %%
        % %
        % % figure(2)
        % % commandprint = sprintf('print -dpdf %salltrajectories',fileonly);
        % % eval(commandprint)
        % %%
        % figure(3)
        % subplot(311)
        % semilogx(frqsFFT_Hz,10*log10(SUUmedian))
        % set(gca,'xlim',[0 6])
        % grid on
        %
        % subplot(312)
        % semilogx(frqsFFT_Hz,10*log10(SRRmedian))
        % set(gca,'xlim',[0 6])
        % grid on
        %
        % subplot(313)
        % semilogx(frqsFFT_Hz,10*log10(abs(SURmedian)))
        % set(gca,'xlim',[0 6])
        % grid on
        pause
    end
end
