%========================== displaySensorRatioSup.m =========================
% this program reads the RatioSups SUT/SREF estimated by spectral approach
% and stored in a drectory as AAresults. That consists on 8 files
% for the 8 sensors of IS26. Each file consists the estimate RatioSups
% on several days of records.
%
% this program calls the program theoreticalSUTresponse4IS26.m
%
clear
addpath ZZtoolbox/
addpath ZZtoolbox/00gabrielson
sensor_UT = 'I26DE_BDF_RSP_2015134_MB3';
saveflag = 0;
% close all
for ihc = 3
    for ii=[2]
        switch ii
            case 1
                comload = sprintf('load AAresultswithFB/resultssta26sensor%i',ihc);
                numfig = ihc+1100;
            case 2
                comload = sprintf('load AAresultswithFB2Ratios/resultssta26sensor%i',ihc);
                numfig = ihc+200;
            case 3
                numfig = ihc+300;
                comload = sprintf('load AAresultswithFBquint/resultssta26sensor%i',ihc);
            case 4
                comload = sprintf('load AAresultswithFBsix/resultssta26sensor%i',ihc);
                numfig = ihc+400;
        end
        
        switch ihc
            case 1
                %                         coeffsens=1.048;
                coeffsens=1.055;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 2
                coeffsens=1.1;%.1;
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
                coeffsens=0.97;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB3';
            case 7
                coeffsens=0.97;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB3';
            case 8
                coeffsens=0.99;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB3';
        end
        
        eval(comload);
        Dstart = str2double(fileswithdotmat(1).name(13:15));
        Dend   = str2double(fileswithdotmat(length(fileswithdotmat)).name(13:15));
        if 1
            doubledaynumber = (Dend-Dstart+3)/2;
        else
            doubledaynumber = 10;
            permutenbmats = randperm(nbmats);
            allRatioSupPfilters = allRatioSupPfilters(:,permutenbmats(1:doubledaynumber));
        end
        
        allRatioSupPfilters_ave          = nanmean(allRatioSupPfilters,2);
        absallRatioSupPfilters_ave       = abs(allRatioSupPfilters_ave);
        angleallRatioSupPfilters_ave_deg = angle(allRatioSupPfilters_ave)*180/pi;
        
        if exist('allRatioInfPfilters','var')
            allRatioInfPfilters_ave          = nanmean(allRatioInfPfilters,2);
            absallRatioInfPfilters_ave       = abs(allRatioInfPfilters_ave);
            angleallRatioInfPfilters_ave_deg = angle(allRatioInfPfilters_ave)*180/pi;
        end
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
        allfrqsPfiltersUZ        = allfrqsPfiltersU(not(allfrqsPfiltersU==0));
        allmeanMSCcstPfiltersU   = allmeanMSCcstPfilters(inda,:);
        nbofvaluesoverthresholdU = nbofvaluesoverthreshold(inda,:);
        allmeanMSCcstPfiltersUZ      = allmeanMSCcstPfiltersU(not(allfrqsPfiltersU==0),:);
        nbofvaluesoverthresholdUZ    = nbofvaluesoverthresholdU(not(allfrqsPfiltersU==0),:);

        allRatioSupPfiltersU         = allRatioSupPfilters(inda,:);
        allRatioSupPfiltersUZ        = allRatioSupPfiltersU(not(allfrqsPfiltersU==0),:);      
        meanallRatioSupPfiltersUZ    = nanmean(allRatioSupPfiltersUZ,2);        
        stdallRatioSupPfiltersUZ     = nanstd(allRatioSupPfiltersUZ,[],2);        
        ICallRatioSupPfiltersUZ      = stdallRatioSupPfiltersUZ ./ ...
            sqrt(sum(nbofvaluesoverthresholdUZ,2));
        
        if exist('allRatioInfPfilters','var')
            allRatioInfPfiltersU         = allRatioInfPfilters(inda,:);
            allRatioInfPfiltersUZ        = allRatioInfPfiltersU(not(allfrqsPfiltersU==0),:);
            meanallRatioInfPfiltersUZ    = nanmean(allRatioInfPfiltersUZ,2);            
            stdallRatioInfPfiltersUZ     = nanstd(allRatioInfPfiltersUZ,[],2);           
            ICallRatioInfPfiltersUZ      = stdallRatioInfPfiltersUZ ./ ...
                sqrt(sum(nbofvaluesoverthresholdUZ,2));
        end

        
        N_freq_vector = 300;
        freq_vector = logspace(log10(0.001),log10(30),N_freq_vector) .';
        [p_total_NRS_sensor, p_total_NRS, TF_ref_sensor, TFsensor4freqRatioSup] = ...
            HCP_acoustical(freq_vector, allfrqsPfiltersUZ, sensor_UT, ref_sensor, 'nofir');
        
        absestimSup = coeffsens * abs(meanallRatioSupPfiltersUZ) .* abs(TFsensor4freqRatioSup);
        absestimInf = coeffsens * abs(meanallRatioInfPfiltersUZ) .* abs(TFsensor4freqRatioSup);
        
        
        %     [FreqFitabs_Hz,absRfit] = smoothpolyLL(allfrqsPfiltersUZ, ...
        %         coeffsens * abs(meanallRatioSupPfiltersUZ),...
        %         segmentsnumber,polydegrees,logtrain,logfit);
        
        figure(numfig)
        clf
        %================================================
        subplot(211)
        semilogx(allfrqsPfiltersUZ,20*log10(absestimSup),'k','linew',1.5),
        
        %== must be very close to the Sup estimates
        if exist('allRatioInfPfilters','var')
            semilogx(allfrqsPfiltersUZ,20*log10(mean([absestimSup, absestimInf],2)),'k','linew',1.5)
            %                     semilogx(allfrqsPfiltersUZ,20*log10(absestimInf),'b','linew',1.5),
        end
        %     hold on
        %         semilogx(ones(2,1)*allfrqsPfiltersUZ',...
        %             [absestim'-ICallRatioSupPfiltersUZ';...
        %             absestim'+ICallRatioSupPfiltersUZ'],'.-b')
        
        %     boxplot(absestim','position',allfrqsPfiltersUZ,'symbol','','whisker',0);
        set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
            'xticklabel',[0.001 0.01 0.1 1 10])
        set(gca,'fontname','times','fontsize',14)
        %     set(gca,'xlim',[0.008 10])
        set(gca,'xlim',[0.015 4.5])
        set(gca,'ylim',2*[-1 1])
        %         set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
        %             'xticklabel',[0.001 0.01 0.1 1 10])
        set(gca,'fontname','times','fontsize',14)
        %     set(gca,'xlim',[0.008 10])
        set(gca,'xlim',[0.01 8])
        set(gca,'ylim',5*[-1 1])
        grid on
        %     xlabel('frequency [Hz]')
        ylabel('Amplitude [dB]','fontname','times','fontsize',14)
        hold on
        plot([1 1]*0.0205,[-45 45],'r','linew',2)
        plot([1 1]*3.98,[-45 45],'r','linew',2)
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
        
        set(gca,'position',[0.1300    0.5056    0.7750    0.3559])
        set(gca,'xtickLabel',[])
        
        ht = text(0.022,-4,'0.02 Hz');
        set(ht,'color','r','fontsize',14,'fontname','times')
        ht = text(2.1,-4,'4 Hz');
        set(ht,'color','r','fontsize',14,'fontname','times')
        ht = text(0.14, -4,'IMS passband');
        set(ht,'color','r','fontsize',14,'fontname','times')
        
        %========================== PHASE =========
        anglestimeSup_rd = angle(meanallRatioSupPfiltersUZ);
        anglestimeInf_rd = angle(meanallRatioInfPfiltersUZ);
        
        anglestimeSup_rd = anglestimeSup_rd + angle(TFsensor4freqRatioSup);
        anglestimeInf_rd = anglestimeInf_rd + angle(TFsensor4freqRatioSup);
        angltheo_rd      = angle(p_total_NRS)+ angle(TF_ref_sensor);
        
        subplot(212)
        
        semilogx(allfrqsPfiltersUZ,unwrap(anglestimeSup_rd)*180/pi,'.-k'),
        if exist('allRatioInfPfilters','var')
            semilogx(allfrqsPfiltersUZ,unwrap(angle(mean([anglestimeSup_rd, anglestimeInf_rd],2)))*180/pi,'k','linew',1.5)
            %                     semilogx(allfrqsPfiltersUZ,20*log10(absestimInf),'b','linew',1.5),
        end
        
        %     boxplot(anglestime_deg','position',allfrqsPfiltersUZ,'symbol','','whisker',0);
        set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
            'xticklabel',[0.001 0.01 0.1 1 10])
        set(gca,'fontname','times','fontsize',14)
        %     set(gca,'xlim',[0.008 10])
        set(gca,'xlim',[0.015 4.5])
        
        set(gca,'ylim',[-20 20])
        
        semilogx(allfrqsPfiltersUZ,unwrap(anglestimeSup_rd)*180/pi,'k','linew',1.5),
        
        %     boxplot(anglestime_deg','position',allfrqsPfiltersUZ,'symbol','','whisker',0);
        %         set(gca,'xscale','log','xtick',[0.001 0.01 0.1 1 10],...
        %             'xticklabel',[0.001 0.01 0.1 1 10])
        set(gca,'fontname','times','fontsize',14)
        %     set(gca,'xlim',[0.008 10])
        set(gca,'xlim',[0.01 8])
        
        set(gca,'ylim',[-40 40])
        
        grid on
        xlabel('frequency [Hz]')
        ylabel('Phase [deg]')
        
        hold on
        
        plot([1 1]*0.0205,[-45 45],'r','linew',2)
        plot([1 1]*3.98,[-45 45],'r','linew',2)
        hold off
        %=============== dipslay
        hold on
        %     semilogx(freq_vector, -unwrap(angltheo_rd)*180/pi, 'r');
        semilogx(freq_vector, unwrap(angltheo_rd)*180/pi-5, 'r--','linew',1.5);
        semilogx(freq_vector, unwrap(angltheo_rd)*180/pi+5, 'r--','linew',1.5);
        hold off
        
        set(gca,'position',[ 0.1300    0.1328    0.7750    0.3333])
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
        fileprintepscmd = sprintf('print -depsc -loose %s3monthsonIS26SUTboxplot%i.eps',printdirectory,ihc);
        fileeps2pdfcmd  = sprintf('!epstopdf %s3monthsonIS26SUTboxplot%i.eps',printdirectory,ihc);
        filermcmd       = sprintf('!rm %s3monthsonIS26SUTboxplot%i.eps',printdirectory,ihc);
        
        eval(fileprintepscmd)
        eval(fileeps2pdfcmd)
        eval(filermcmd)
        if saveflag
            fileprintpngcmd = sprintf('print -dpng -loose ../../textes/6distConjointHMSC/figures/3monthsonIS26SUT%i',ihc);
            fileprintpngcmd = sprintf('print -dpdf -loose threemonthsonIS26SUT%i',ihc);
            eval(fileprintpngcmd)
        end
        
        
    end
    
    if saveflag
        
        %             eval(fileprintepscmd)
        %             eval(fileeps2pdfcmd)
        %             eval(filermcmd)
        
        %                 for alfred
        fileprintjpgcmd = sprintf('print -dpng -loose ../../../../3monthsonIS26SUT%i%i',ihc,ii);
        eval(fileprintjpgcmd)
        
        
    end
end
