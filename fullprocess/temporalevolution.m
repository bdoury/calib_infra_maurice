% ratio at 1 Hz, look for outliers
%=================================
clear
addpath ZZtoolbox/
addpath ZZtoolbox/00gabrielson
close all
saveflag = 1;
bootdraw = 0;
nbdraw=100;
Ncouples=50;
SUFFICIENTNUMBER = 2;
for ihc=8
    % keep 2
    for ii=1
        switch ii
            case 1
                comload = sprintf('load AAresultsaround1Hz/resultssta26sensor%i',ihc);
                numfig = ihc+100;
            case 2
                comload = sprintf('load AAresults0812Hzbis/resultssta26sensor%i',ihc);
                numfig = ihc+200;
            case 3
                numfig = ihc+300;
                comload = sprintf('load AAresults0614Hz/resultssta26sensor%i',ihc);
        end
        switch ihc
            case 1
                coeffsens=1.02;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 2
                coeffsens=1.06;%.1;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 3
                coeffsens=1.032;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 4
                coeffsens=1.07;
                ref_sensor = 'I26DE_BDF_RSP_2015134_MB2005';
            case 5
                coeffsens=1.05;
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
        
        filtercharact.SCPperiod_sec
        
        Dstart = str2double(fileswithdotmat(1).name(13:15));
        Dend   = str2double(fileswithdotmat(length(fileswithdotmat)).name(13:15));
        %         if 1
        %             doubledaynumber = (Dend-Dstart+3)/2;
        %         else
        %             doubledaynumber = 10;
        %             permutenbmats = randperm(nbmats);
        %             allRatioSupPfilters = allRatioSupPfilters(:,permutenbmats(1:doubledaynumber));
        %         end
        
        sensor_UT = 'I26DE_BDF_RSP_2015134_MB3';
        
        N_freq_vector = 1;
        freq_vector   = 1;
        [p_total_NRS_sensor, p_total_NRS, TF_ref_sensor, TFsensor4freqRatio] = ...
            HCP_acoustical(freq_vector, freq_vector, sensor_UT, ref_sensor, 'nofir');
        expectedvalueat1Hz = abs(p_total_NRS_sensor);
        
        indselect = find(allfrqsPfilters(:,1)>=1,1,'first');
        
        allpraticalvaluesat1Hz = coeffsens*abs(TFsensor4freqRatio)*abs(allRatioSupPfilters(indselect,:)) ./ ...
            expectedvalueat1Hz;
        
        nbofcouplesdays = length(allpraticalvaluesat1Hz);
        %===================== plots ==============
        figure(numfig)
        subplot(211);
        
        if bootdraw
            if 0
                for idraw=0:fix(nbofcouplesdays/Ncouples)-1
                    allpraticalvaluesat1Hzrand = allpraticalvaluesat1Hz(:,randperm(nbofcouplesdays));
                    id1=idraw*Ncouples+1;
                    id2=id1+Ncouples-1;
                    meanallpraticalvaluesat1Hz_ii=nanmean(allpraticalvaluesat1Hzrand(:,id1:id2));
                    plot(idraw,20*log10(meanallpraticalvaluesat1Hz_ii),'ok','markersize',6,'markerfacec','k')
                    hold on
                end
            else
                for idraw=1:nbdraw
                    PP_ii = randperm(nbofcouplesdays);
                    meanallpraticalvaluesat1Hz_ii=nanmean(allpraticalvaluesat1Hz(:,PP_ii(1:Ncouples)));
                    
                    plot(idraw,20*log10(meanallpraticalvaluesat1Hz_ii),'ok','markersize',6,'markerfacec','k')
                    
                    hold on
                    %             sigmaonRatio = allSTDmodRatioPfilters(indselect,:) ./ ...
                    %                 sqrt(nbofvaluesoverthreshold(indselect,:));
                    %         plot(ones(2,1)*(1:size(allRatioSupPfilters,2)),...
                    %             coeffsens*[abs(allRatioSupPfilters(indselect,:))-sigmaonRatio;...)
                    %             abs(allRatioSupPfilters(indselect,:))+sigmaonRatio],'.-','color',0.6*[1 1 1]);
                end
            end
        else
            
            indnotsufficient = find(nbofvaluesoverthreshold(indselect,:)<=SUFFICIENTNUMBER);
            allpraticalvaluesat1Hz(indnotsufficient) = NaN;
            plot(20*log10(allpraticalvaluesat1Hz),'ok','markersize',6,'markerfacec','k')
            hold on
        end
        
        
        %         plot([0 10000], [0 0],'--r','linew',2)
        plot([0 10000], 20*log10(1.05*[1 1]),'--r','linew',2)
        plot([0 10000], 20*log10(0.95*[1 1]),'--r','linew',2)
        set(gca,'fontname','times','fontsize',14)
        ylabel('Gain at 1 Hz','fontname','times','fontsize',14)
        grid on
        if bootdraw
            set(gca,'xlim',[0 idraw]);
        else
            set(gca,'xlim',[1 nbofcouplesdays])
        end
        set(gca,'ylim',[-1.5 1.5])
        
        title(sprintf('IS26 -  sensor #%i, threshold = %4.2f\nday number = %i, T = %i s',...
            ihc, MSCthreshold, 2*nbofcouplesdays, filtercharact.SCPperiod_sec),'fontname','times','fontsize',14)
        
        subplot(212);
        %             plot(allmeanMSCcstPfilters(indselect,:),'ob','markersize',6,...
        %                 'markerfacec','b')
        semilogy(nbofvaluesoverthreshold(indselect,:),'ob','markersize',6,...
            'markerfacec','b')
        %             set(gca,'ylim',[0.98 1])
        ylabel('number fo values','fontname','times','fontsize',14)
        grid on
        set(gca,'fontname','times','fontsize',14)
        if bootdraw
            set(gca,'xlim',[0 idraw]);
        else
            set(gca,'xlim',[1 nbofcouplesdays])
        end
        
        
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
        fileprintepscmd = sprintf('print -depsc -loose %sevolutionon%iatfreq1bis.eps',printdirectory,ihc);
        fileeps2pdfcmd  = sprintf('!epstopdf %sevolutionon%iatfreq1.eps',printdirectory,ihc);
        filermcmd       = sprintf('!rm %sevolutionon%iatfreq1.eps',printdirectory,ihc);
        %
        if saveflag
            eval(fileprintepscmd)
            %             eval(fileeps2pdfcmd)
            %             eval(filermcmd)
        end
    end
    
end
