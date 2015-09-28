function coeffsens = theoreticalNRSonIS26(Rfit, Frqs4fit_Hz, RSA, FreqsSA_Hz, ihc, figdisplay)


% used function : 
%        - IS26TrosetteBIS.m
%=========================================================
%
%  Cavity specs (microbarometer, primary, secondary,
%                     and tertiary manifolds)
%=========================================================
%  Primary pipe characteristics
NRS.Npri  = 4;%           --- number
NRS.Nsec  = 24;%          --- number

NRS.Lpri  = 5.75; %10.36;  --- length in meter
NRS.Lsec  = 3.25; %3;%14.97;     --- length in meter

NRS.a_pri = 0.824*25.4e-3/2; % 13e-3/2;%    --- radius in meter
NRS.a_sec = 0.622*25.4e-3/2; % 13e-3/2;%   --- radius in meter

NRS.L_stub = 1.6;%3.708;%
NRS.a_stub = 0.622*25.4e-3/2; %13e-3/2;% 36-m HCP

NRS.a_pri_cavity = 168e-3/2;%50e-3;% --- radius in meter meter
NRS.L_pri_cavity = 40e-3;%38e-3;%    --- height in

NRS.a_sec_cavity = 100e-3/2;%168e-3/2;%    --- radius in meter
NRS.L_sec_cavity = 38e-3; %80e-3;%    --- height in meter

NRS.a_microbarometer = 70e-3;  %  radius and length to get correct MB2000
NRS.L_microbarometer = 39e-3;   %  volume (0.0006 m^3 from Alcoverro)
% Length of stub from primary manifold to microbarometer [m]
% WARNING: the smallest wavelength is related to the high frequency
% typically the 300/10 about 30 meters
%
Ninlet = NRS.Npri * NRS.Nsec;
%=========================================================
% Gabrielson code computes the frequency response
TFNRSfit              = Ninlet*IS26TrosetteBIS(Frqs4fit_Hz, NRS);
%=========================================================

%=========================================================
% Benoit code computes the frequency response of the
% theoretical sensors:
responseMB3theo    = idc2fap_is(FreqsSA_Hz,'I26DE_BDF_RSP_2015134_MB3');
responseMB2005theo = idc2fap_is(FreqsSA_Hz,'I26DE_BDF_RSP_2015134_MB2005');
%=========================================================
responseMB3theo_fit    = idc2fap_is(Frqs4fit_Hz,'I26DE_BDF_RSP_2015134_MB3');
responseMB2005theo_fit = idc2fap_is(Frqs4fit_Hz,'I26DE_BDF_RSP_2015134_MB2005');
%=========================================================

RSA_ave  = nanmean(RSA,2);

nbmats = size(RSA,2);

if ihc<6
    absresponseSUTtheo           = abs(TFNRSfit)  .* abs(responseMB2005theo_fit);
    absresponseSUTtheo_sup       = 1.05* (abs(TFNRSfit)  .* abs(responseMB2005theo_fit));
    absresponseSUTtheo_inf       = 0.95* (abs(TFNRSfit)  .* abs(responseMB2005theo_fit));
    absresponseSUTestim          = abs(Rfit) .* abs(responseMB2005theo_fit);
    absallRatioPfilters_ave_corr = abs(RSA_ave .* responseMB2005theo);
    absallRatioPfilters_corr     = abs(RSA) .* (abs(responseMB2005theo) * ones(1,nbmats));
    
    angleresponseSUTtheo         = 180*(angle(TFNRSfit)  + angle(responseMB2005theo_fit))/pi;
    angleresponseSUTestim        = (angle(Rfit)+ angle(responseMB2005theo_fit))*180/pi;
    angleallRatioPfilters_ave_deg_corr = (angle(RSA_ave)+angle(responseMB2005theo))*180/pi;
    
else
    absresponseSUTtheo           = abs(TFNRSfit) .* abs(responseMB3theo_fit);
    absresponseSUTtheo_sup       = 1.05* (abs(TFNRSfit) .* abs(responseMB3theo_fit));
    absresponseSUTtheo_inf       = 0.95* (abs(TFNRSfit) .* abs(responseMB3theo_fit));
    absresponseSUTestim          = abs(Rfit) .* abs(responseMB3theo_fit);
    absallRatioPfilters_ave_corr = abs(RSA_ave .* responseMB3theo);
    absallRatioPfilters_corr     = abs(RSA) .* (abs(responseMB3theo) * ones(1,nbmats));
    
    angleresponseSUTtheo         = 180*(angle(TFNRSfit)  + angle(responseMB3theo_fit))/pi;
    angleresponseSUTestim        = (angle(Rfit)+ angle(responseMB3theo_fit))*180/pi;
    angleallRatioPfilters_ave_deg_corr = (angle(RSA_ave)+angle(responseMB3theo))*180/pi;
    
end


figure(figdisplay)

switch ihc
    case 1
        coeffsens = 1.05;
    case 2
        coeffsens = 1.1;
    case 3
        coeffsens = 1.07;
    case 4
        coeffsens = 1.07;
    case 5
        coeffsens = 1.1;
    case 6
        coeffsens = 1.01;
    case 7
        coeffsens = 0.97;
    case 8
        coeffsens = 1.02;
end

subplot(222)
semilogx(Frqs4fit_Hz, (absresponseSUTtheo),'r')
hold on
semilogx(Frqs4fit_Hz, (absresponseSUTtheo_sup),'r--')
semilogx(Frqs4fit_Hz, (absresponseSUTtheo_inf),'r--')
semilogx(Frqs4fit_Hz, (absresponseSUTestim*coeffsens),'m') ; %0.94, 1.02
grid on
hold off

set(gca,'xlim',[0.002 8])
set(gca,'ylim',[0.9 1.1])
set(gca,'fontname','times','fontsize',10)

subplot(224)
semilogx(Frqs4fit_Hz, (angleresponseSUTtheo),'r')
hold on
semilogx(Frqs4fit_Hz, (angleresponseSUTestim),'m') ; %0.94, 1.02
grid on
hold off
set(gca,'fontname','times','fontsize',10)

set(gca,'xlim',[0.002 8])
set(gca,'ylim',[-10 10])

subplot(221)
semilogx(FreqsSA_Hz,coeffsens*absallRatioPfilters_corr,'.', ...
    'color',0.9*ones(3,1))
hold on
semilogx(FreqsSA_Hz,coeffsens*(absallRatioPfilters_ave_corr),'.')
hold on
semilogx(Frqs4fit_Hz, (absresponseSUTestim*coeffsens),'m','linew',2)
hold off
set(gca,'fontname','times','fontsize',10)
set(gca,'xlim',[0.002 8])
set(gca,'ylim',[0.9 1.1])
grid on
xlabel('frequency - Hz')
ylabel('Estimated Gain')
drawnow
%=====
subplot(223)
%     semilogx(allfrqsPfilters,angle(allRatioPfilters)*180/pi,'.', ...
%         'color',0.9*ones(3,1))
%     hold on
semilogx(FreqsSA_Hz,angleallRatioPfilters_ave_deg_corr,'.')
hold on
semilogx(Frqs4fit_Hz, (angleresponseSUTestim),'m') ; %0.94, 1.02
hold off
set(gca,'fontname','times','fontsize',10)
set(gca,'xlim',[0.002 8])
set(gca,'ylim',[-10 10])
grid on
xlabel('frequency - Hz')
ylabel('Estimated Phase - degree')
drawnow
%=============================================================================