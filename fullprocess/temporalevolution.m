clear
addpath ZZtoolbox/
addpath ZZtoolbox/00gabrielson
close all
saveflag = 0;
for ihc=1
    % keep 2
    for ii=2
        switch ii
            case 1
                comload = sprintf('load AAresults0812Hz/resultssta26sensor%i',ihc);
                numfig = ihc+100;
            case 2
                comload = sprintf('load AAresults0812Hzbis/resultssta26sensor%i',ihc);
                numfig = ihc+200;
            case 3
                numfig = ihc+300;
                comload = sprintf('load AAresults0812Hzter/resultssta26sensor%i',ihc);
        end
        switch ihc
            case 1
                coeffsens=1.035;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 2
                coeffsens=1.06;%.1;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 3
                coeffsens=1.052;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 4
                coeffsens=1.06;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 5
                coeffsens=1.065;
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
        
        sensor_UT = 'I26DE_BDF_RSP_2015134_MB3';
        
        
        N_freq_vector = 1;
        freq_vector   = 1;
        [p_total_NRS_sensor, p_total_NRS, TF_ref_sensor, TFsensor4freqRatio] = ...
            HCP_acoustical(freq_vector, freq_vector, sensor_UT, ref_sensor, 'nofir');
        expectedvalueat1Hz = abs(p_total_NRS_sensor);
        
        indselect = find(allfrqsPfilters(:,1)>=1,1,'first');
        
        allpraticalvaluesat1Hz = coeffsens*abs(TFsensor4freqRatio)*abs(allRatioPfilters(indselect,:)) ./ ...
            expectedvalueat1Hz;
        
        %===================== plots ==============
        figure(numfig)
        subplot(211);
        plot(allpraticalvaluesat1Hz,'ok','markersize',6,'markerfacec','k')
        
        hold on
        sigmaonRatio = allSTDmodRatioPfilters(indselect,:) ./ ...
            sqrt(nbofvaluesoverthreshold(indselect,:));
        %         plot(ones(2,1)*(1:size(allRatioPfilters,2)),...
        %             coeffsens*[abs(allRatioPfilters(indselect,:))-sigmaonRatio;...)
        %             abs(allRatioPfilters(indselect,:))+sigmaonRatio],'.-','color',0.6*[1 1 1]);
        
        
        
        plot([0 100], [1 1],'--r','linew',2)
        plot([0 100], 1.05*[1 1],'--r','linew',2)
        plot([0 100], 0.95*[1 1],'--r','linew',2)
        hold off
        set(gca,'fontname','times','fontsize',14)
        ylabel('Gain at 1 Hz','fontname','times','fontsize',14)
        grid on
        set(gca,'xlim',[1 ceil(doubledaynumber)])
        set(gca,'ylim',[0.9 1.12])
        
        title(sprintf('IS26 -  sensor #%i, threshold = %4.2f\nday number = %i, T = %i s',...
            ihc, MSCthreshold, 2*doubledaynumber, filtercharact.SCPperiod_sec),'fontname','times','fontsize',14)
        
        subplot(212);
        %             plot(allmeanMSCcstPfilters(indselect,:),'ob','markersize',6,...
        %                 'markerfacec','b')
        semilogy(nbofvaluesoverthreshold(indselect,:),'ob','markersize',6,...
            'markerfacec','b')
        %             set(gca,'ylim',[0.98 1])
        ylabel('number fo values','fontname','times','fontsize',14)
        grid on
        set(gca,'fontname','times','fontsize',14)
        
        
        HorizontalSize = 16;
        VerticalSize   = 11;
        set(gcf,'units','centimeters');
        set(gcf,'paperunits','centimeters');
        set(gcf,'PaperType','a3');
        set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
        set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
        set(gcf,'color', [1,1,0.92]);
        set(gcf, 'InvertHardCopy', 'off');
        
        printdirectory  = ' ../slidesITW2015/';
        fileprintepscmd = sprintf('print -depsc -loose %sevolutionon%iatfreq1.eps',printdirectory,ihc);
        fileeps2pdfcmd  = sprintf('!epstopdf %sevolutionon%iatfreq1.eps',printdirectory,ihc);
        filermcmd       = sprintf('!rm %sevolutionon%iatfreq1.eps',printdirectory,ihc);
        %
        if saveflag
            eval(fileprintepscmd)
            eval(fileeps2pdfcmd)
            eval(filermcmd)
        end
    end
end
