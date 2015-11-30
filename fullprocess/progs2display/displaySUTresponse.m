%========================== displaySUTresponse.m =========================
% this program reads the ratios SUT/SREF estimated by spectral approach
% and stored in a drectory as AAresults. That consists on 8 files
% for the 8 sensors of IS26. Each file consists the estimate ratios
% on several days of records.
%=========
% this program calls the function:
%         - HCP_acoustical.m
% using both ZZtoolbox/00gabrielson and ../ZZtoolbox/00benoit
%
%=========================================================================
clear
addpath ../ZZtoolbox/
addpath ../ZZtoolbox/00gabrielson
addpath ../ZZtoolbox/00benoit
directorysignals    = '../../../../AAdataI26calib/';

%==== this directory contains the parameters evalauted by the
% program estimationwithFB.m
directoryinputresults = '../AAresultswithFB98/';

sensor_UT    = 'I26DE_BDF_RSP_2015134_MB3';
saveflag     = 1;
trimmeanflag = 0;
remainindex1 = [1:61  63:70];

for ihc = 1
    % list of the files from 1 to nbmats
    % if you want a name type fileswithdotmat(#)
    fileswithdotmat = dir(sprintf('../%ss%i/s%iy*.mat',...
        directorysignals,ihc,ihc));
    comload = sprintf('load %sresultssta26sensor%i',directoryinputresults,ihc);
    eval(comload);
    switch ihc
        case 1
            remainindex = remainindex1; % 2015/10/13
        case 2
            remainindex = [1:nbmats]; % 2015/10/13
        case 4
            remainindex = [1:65 66 67:nbmats]; % 2015/10/11
        otherwise
            remainindex = [1:nbmats];
    end
    switch ihc
        case 1
            coeffsens=1.04;
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
            coeffsens=1;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB3';
        case 7
            coeffsens=0.97;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB3';
        case 8
            coeffsens=1.02;
            ref_sensor = 'I26DE_BDF_RSP_2015134_MB3';
    end
    %%
    nbdoubledays                 = length(remainindex);
    allRatioSupPfilters          = allRatioSupPfilters(:,remainindex);
    allSTDmodRatioSupPfilters    = allSTDmodRatioSupPfilters(:,remainindex);
    allSTDphaseRatioSupPfilters  = allSTDphaseRatioSupPfilters(:,remainindex);
    
    allRatioInfPfilters          = allRatioInfPfilters(:,remainindex);
    allSTDmodRatioInfPfilters    = allSTDmodRatioInfPfilters(:,remainindex);
    allSTDphaseRatioInfPfilters  = allSTDphaseRatioInfPfilters(:,remainindex);
    allmeanMSCcstPfilters        = allmeanMSCcstPfilters(:,remainindex);
    nbofvaluesoverthreshold      = nbofvaluesoverthreshold(:,remainindex);
    allScpPfilters               = allScpPfilters(:,:,remainindex);


    STDmodRatioPfilters_ave        = nanmedian(allSTDmodRatioSupPfilters,2);
    STDphaseRatioPfilters_ave      = nanmedian(allSTDphaseRatioSupPfilters,2);
    
    %== sort in increasing order
    [allfrqsPfiltersS, inds]       = sort(allfrqsPfilters);
    allRatioPfiltersS              = allRatioSupPfilters(inds,:);
    allmeanMSCcstPfiltersS         = allmeanMSCcstPfilters(inds,:);
    STDmodRatioPfilters_aveS       = STDmodRatioPfilters_ave(inds);
    STDphaseRatioPfilters_aveS     = STDphaseRatioPfilters_ave(inds);
    nbofvaluesoverthresholdS       = nbofvaluesoverthreshold(inds,:);
    
    %== unique frequency value
    [allfrqsPfiltersUS, inda]      = unique(allfrqsPfiltersS);
    allRatioPfiltersUS             = allRatioPfiltersS(inda,:);
    allmeanMSCcstPfiltersUS        = allmeanMSCcstPfiltersS(inda,:);
    STDmodRatioPfilters_aveUS      = STDmodRatioPfilters_aveS(inda);
    STDphaseRatioPfilters_aveUS    = STDphaseRatioPfilters_aveS(inda);
    nbofvaluesoverthresholdUS      = nbofvaluesoverthresholdS(inda,:);
    
    %== without 0 frequency values
    allfrqsPfiltersUSZ             = allfrqsPfiltersUS(not(allfrqsPfiltersUS==0));
    RatioPfiltersUSZ               = allRatioPfiltersUS(not(allfrqsPfiltersUS==0),:);
    allmeanMSCcstPfiltersUSZ       = allmeanMSCcstPfiltersUS(not(allfrqsPfiltersUS==0),:);
    STDmodRatioPfilters_aveUSZ     = STDmodRatioPfilters_aveUS(not(allfrqsPfiltersUS==0));
    STDphaseRatioPfilters_aveUSZ   = STDphaseRatioPfilters_aveUS(not(allfrqsPfiltersUS==0));
    nbofvaluesoverthresholdUSZ     = nbofvaluesoverthresholdUS(not(allfrqsPfiltersUS==0),:);
    
    %====== absolute and arg of the ratios
    modRatioPfiltersUSZ            = abs(RatioPfiltersUSZ);
    phaseRatioPfiltersUSZ_rd       = angle(RatioPfiltersUSZ);
    %====== averaging by TRIMMEAN to avoid outliers
    if trimmeanflag
        meanmodRatioPfiltersUSZ        = trimmean(modRatioPfiltersUSZ,30,2);
        meanphasePfiltersUSZ_rd        = trimmean(phaseRatioPfiltersUSZ_rd,30,2);
    else
        meanmodRatioPfiltersUSZ        = nanmean(modRatioPfiltersUSZ,2);
        meanphasePfiltersUSZ_rd        = nanmean(phaseRatioPfiltersUSZ_rd,2);
    end
    %====== STDs
    STDmodPfiltersUSZ              = nanstd(modRatioPfiltersUSZ,[],2);
    STDphasePfiltersUSZ_rd         = nanstd(phaseRatioPfiltersUSZ_rd,[],2);
    
    ICallRatioPfiltersUSZ          = STDmodPfiltersUSZ ./ ...
        sqrt(sum(nbofvaluesoverthresholdUSZ,2));
    
    ICallRatioPfiltersUSZbis       = STDmodRatioPfilters_aveUSZ ./ ...
        sqrt(sum(nbofvaluesoverthresholdUSZ,2));
    
    N_freq_vector = 300;
    freq_vector   = logspace(log10(0.001),log10(30),N_freq_vector) .';
    [p_total_NRS_sensor, p_total_NRS, TF_ref_sensor, TFsensor4freqRatio] = ...
        HCP_acoustical(freq_vector, allfrqsPfiltersUSZ, sensor_UT, ref_sensor, 'nofir');
    
    %================================
    absestimwithcorrect = coeffsens * meanmodRatioPfiltersUSZ .* abs(TFsensor4freqRatio);
    %================================================
    subplot(211)
    semilogx(allfrqsPfiltersUSZ,20*log10(absestimwithcorrect),'ko-',...
        'markersize',4,'markerfaceco','k'),
    grid on
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
    title(sprintf('IS26 -  sensor H%i, MSC threshold = %4.2f\nday number = %i',...
        ihc, MSCthreshold, 2*nbdoubledays),'fontname','times','fontsize',14)
    %         title(sprintf('IS26 -  sensor H%i\ndashed line: +/-5%s for amplitude, +/- 5 degrees for phase', ihc,'%'),...
    %             'fontname','times','fontsize',14)
    
    set(gca,'fontname','times','fontsize',14)
    set(gca,'xlim',[0.01 5])
    set(gca,'ylim',4*[-1 1])
    set(gca,'xtickLabel',[])
    %              xlabel('frequency [Hz]')
    
    set(gca,'position',[0.1300    0.5056    0.7750    0.3559])
    
    ht = text(0.022,-3.4,'0.02 Hz');
    set(ht,'color','r','fontsize',14,'fontname','times')
    ht = text(2.1,-3.4,'4 Hz');
    set(ht,'color','r','fontsize',14,'fontname','times')
    ht = text(0.14, -3.4,'IMS passband');
    set(ht,'color','r','fontsize',14,'fontname','times')
    
    %========================== PHASE =========
    
    anglewithcorrect_rd = meanphasePfiltersUSZ_rd + angle(TFsensor4freqRatio);
    angltheo_rd         = angle(p_total_NRS)+ angle(TF_ref_sensor);
    
    subplot(212)
    semilogx(allfrqsPfiltersUSZ,unwrap(anglewithcorrect_rd)*180/pi,...
        'ko-','markersize',4,'markerfaceco','k')
    
    set(gca,'fontname','times','fontsize',14)
    
    grid on
    xlabel('frequdoubledaynumberency [Hz]')
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
    
    set(gca,'fontname','times','fontsize',14)
    set(gca,'xlim',[0.01 5])
    set(gca,'ylim',40*[-1 1])
    xlabel('frequency [Hz]')
    
    set(gca,'position',[ 0.1300    0.1328    0.7750    0.3333])
    %==============================================================
    HorizontalSize = 16;
    VerticalSize   = 14;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a3');
            set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    set(gcf,'color', [1,1,0.92]);
    set(gcf, 'InvertHardCopy', 'off');
    
    printdirectory  = ' ../../figures/';
    
    
    if trimmeanflag==1
        fileprint = sprintf('%swithproblemons%iSolved.eps',printdirectory,ihc);
    elseif nbdoubledays<nbmats
        fileprint = sprintf('%swithoutproblemons%i.eps',printdirectory,ihc);
    else
        fileprint = sprintf('%swithproblemons%i.eps',printdirectory,ihc);
    end
    
    fileprintepscmd = sprintf('print -depsc -loose %s',fileprint);
    fileeps2pdfcmd  = sprintf('!epstopdf %s',fileprint);
    filermcmd       = sprintf('!rm %s',fileprint);
    
end
figure(1)
if saveflag
    eval(fileprintepscmd)
%     eval(fileeps2pdfcmd)
%     eval(filermcmd)
end
%%

% ip=1;
% NaverageFFTs = filtercharact(1).ratioDFT2SCP;
% allT.TUUonUR = linspace(0.7,1.3,100);
% allT.TURonRR = linspace(0.7,1.3,100);
% allT.MSC     = linspace(0.5,1,100);
% allT.phase   = linspace(0,2*pi,100);
%
% listifq = cumsumnbfq_ip(ip,1):cumsumnbfq_ip(ip,2);
% Llistifq = length(listifq);
% STDmodtheo_ip = zeros(Llistifq,1);
%
% for indfq=1:Llistifq
%     ifq=listifq(indfq);
%     SCP_ip_ifq = nanmean(allScpPfilters(:,ifq,:),3);
%     if any(isnan(SCP_ip_ifq))
%         STDmodtheo_ip(indfq)=NaN;
%     else
%         RR=[SCP_ip_ifq(1) SCP_ip_ifq(3);SCP_ip_ifq(3)' SCP_ip_ifq(2)];
%         [statUUonUR, statURonRR, statMSC]    = ...
%             statsRatiosHbis(allT, RR, NaverageFFTs, 0.3);
%         STDmodtheo_ip(indfq) = diff(statUUonUR.CI)/2;
%     end
% end
% nbofvalues_ip = sum(nbofvaluesoverthreshold(cumsumnbfq_ip(ip,1):cumsumnbfq_ip(ip,2),:),2);
% STDmodempiric_ip = nanmean(...
%     allSTDmodRatioSupPfilters(cumsumnbfq_ip(ip,1):cumsumnbfq_ip(ip,2),:),2);
% corrlevel = corr(STDmodtheo_ip(and(not(isnan(STDmodtheo_ip)),not(STDmodempiric_ip==0))), STDmodempiric_ip(and(not(isnan(STDmodtheo_ip)),not(STDmodempiric_ip==0))));
%
% [STDmodtheo_ip./(STDmodempiric_ip) nbofvalues_ip]
% corrlevel
%
% figure(1)
% subplot(121)
% plot( allT.TUUonUR,statURonRR.pdf,'.-')
% hold on
% plot( allT.TUUonUR,statUUonUR.pdf,'.-r')
% hold off
% subplot(122)
% plot(STDmodtheo_ip, STDmodempiric_ip,'.')
%
%
