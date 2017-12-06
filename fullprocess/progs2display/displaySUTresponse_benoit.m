%========================== displaySUTresponse.m =========================
% this program displays the NRS response at IS26
%=========
% this program calls the function:
%         - HCP_acoustical_benoit.m
% using both ZZtoolbox/00gabrielson and ../ZZtoolbox/00benoit
%
%=========================================================================
clear TF_ref_sensor TFsensor4freqRatio angltheo_rd angltheo_rd_NRS directoryinputresults ihc numfig ref_sensor saveprintflag sensor_UT trimmeanflag;

clear N_freq_vector freq_vector p_total_NRS p_total_NRS_sensor;
addpath ../ZZtoolbox/
addpath ../ZZtoolbox/00gabrielson
addpath ../ZZtoolbox/00benoit

%==== this directory contains the parameters evalauted by the
% program estimationwithFB.m

sensor_UT     = 'I26DE_BDF_RSP_2015134_MB3';

for ihc = 1:1
    numfig = ihc+1;
    figure(numfig);
    
    N_freq_vector = 999;
    freq_vector   = logspace(-7,2,N_freq_vector) .';
    [p_total_NRS_sensor, p_total_NRS] = ...git status
        HCP_acoustical_benoit(freq_vector, sensor_UT, 'nofir');

    subplot(211)

    semilogx(freq_vector, 20*log10(abs(p_total_NRS)), 'r');
    hold on;
    %semilogx(freq_vector, 20*log10(abs(p_total_NRS_sensor)), 'g');
    grid on
    ylabel('Amplitude [dB]','fontname','times','fontsize',14)
    hold off

    set(gca,'fontname','times','fontsize',14)
    set(gca,'xlim',[10E-7 60])
    set(gca,'ylim',25*[-0.5 1])
    set(gca,'xtickLabel',[])
    
    set(gca,'position',[0.1300    0.5056    0.7750    0.3559])

    
    %========================== PHASE =========
   
    angltheo_rd         = angle(p_total_NRS_sensor);
    angltheo_rd_NRS         = angle(p_total_NRS);
        %=============== dipslay
    subplot(212)
    semilogx(freq_vector, unwrap(angltheo_rd_NRS)*180/pi, 'r');
    hold on
    %semilogx(freq_vector, unwrap(angltheo_rd)*180/pi, 'g');
    set(gca,'fontname','times','fontsize',14)
    
    grid on
    xlabel('frequdoubledaynumberency [Hz]')
    ylabel('Phase [deg]')
    hold off


    set(gca,'fontname','times','fontsize',14)
    set(gca,'xlim',[10E-7 60])
    set(gca,'ylim',80*[-1 1])
    xlabel('frequency [Hz]')
    
end
