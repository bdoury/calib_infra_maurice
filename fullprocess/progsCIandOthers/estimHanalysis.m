%===============================================================
% estimHanalysis.m
%===============================================================
% Determine the level of MSC useful to reach a given level of
% accuracy for the ratio HUT/HREF.
% Simulation is conducted on synthetical or semi-synthetical
% signals.
%===============================================================
%
% used function: fbankanalysis.m
% with no filterbank
%
% We draw randomly and use the H estimate
% Models of signals
% xREF = BGfield + sqrt(sigma2REF)*wREF;
% xUT  = BGfield + sqrt(sigma2UT)*wUT;
% yREF = filter(bREF,aREF,xREF);
% yUT  = filter(bUT,aUT,xUT);
%
% Results have to be compared to these
% of program STATSONH11H12.m, plotted with allprintstatsonHest.m
%===============================================================
% N is divided in segments of length equal to Tfft_sec*Fs_Hz
% with Tfft_sec=300.
% To form the data block on which the spectrum is performed, 
% we use M consecutive segments. 
% The simulation is iterated in such a way the duration of the
% total observation is 5 hours. This value is obvioulsy not crucial.
% It is just to perform the RMSE on enough simulations.
% The RMSE is computed by averaging on the full frequency band.
%
%===============================================================
clear
addpath  ../../../ZZtoolbox/
synth            = 1;
Fs_Hz   = 20;
FFT_sec = 300;
Lfft    = fix(FFT_sec*Fs_Hz);
duration_hour    = 5;
N           = duration_hour*3600*Fs_Hz;
frqs           = (0:Lfft-1)*Fs_Hz/Lfft;
selectbandave_Hz = 10;
idbband          = find(frqs<selectbandave_Hz,1,'last');
bandave          = 1:idbband;


%======== filter REF ==========
FREF_Hz      = 0;
rhoREF       = 0.7;
costhetaREF  = (1+rhoREF*rhoREF)*cos(2*pi*FREF_Hz/Fs_Hz)/2/rhoREF;
bREF         = 0.5*[1+rhoREF*rhoREF;-4*rhoREF*costhetaREF;1+rhoREF*rhoREF];
aREF         = [1 -2*rhoREF*costhetaREF rhoREF*rhoREF]';
HREF         = fft(bREF,Lfft) ./ fft(aREF,Lfft);
HREF         = HREF(:);
%======== filter UT ==========
FUT_Hz       = 0;
rhoUT        = 0.7;
costhetaUT   = (1+rhoUT*rhoUT)*cos(2*pi*FUT_Hz/Fs_Hz)/2/rhoUT;
bUT          = 0.5*[1+rhoUT*rhoUT; -4*rhoUT*costhetaUT; 1+rhoUT*rhoUT];
aUT          = [1 -2*rhoUT*costhetaUT rhoUT*rhoUT]';
HUT          = fft(bUT,Lfft) ./ fft(aUT,Lfft);
HUT          = HUT(:);

% randn('seed',0)
if synth
    % for synthetical signals
    % the BGfield is white
    randnN = randn(N,1);
    randnN = randnN/std(randnN);
else
    % the BGfield is drawn
    % from the database
    directorydata = '../../../../DATA_IS/I26/';
    % directorydata = './';
    filesmat = dir([directorydata '*.mat']);
    nbmats   = length(filesmat);
    ifile    = 2;
    commandload = sprintf('load %s%s',directorydata,filesmat(ifile).name);
    eval(commandload)
    id1=fix(rand*1e6);
    data_pa = records{2}.data(id1+(1:N));
    data_pa = data_pa-mean(data_pa);
    randnN = data_pa/std(data_pa);
end
% theoretically gamma has no effect
gamma          = 5;
BGfield        = sqrt(gamma)*randnN;

overlapFFT     = 0.5;
overlapAVE     = 0;

sqrtLfft       = sqrt(Lfft);
frqs           = (0:Lfft-1)*Fs_Hz/Lfft;

listM          = [5 10];
LM             = length(listM);
listMSC        = 0.8:0.02:0.99;
LMSC           = length(listMSC);
% root mean square error on the ratio SUU/abs(SUR)
RMSEUUUR       = zeros(LMSC,LM);
RMSEAB2        = zeros(LMSC,LM);

%==============================
% Because airinlet is given the MSC implies
% the level of white noise (not true for non white)
%==============================
for iMSC=1:LMSC
    MSCtrue   = listMSC(iMSC);
    %==============================
    % airinlet induces lower noise on UT than on REF
    % airinlet value is of order of magnitude of the
    % number of airinlets, typically a few tens
    % but has little consequence because we fix the MSC
    % by taking into  account this number.
    %==============================
    airinlet    = 48;
    sigma2UT    = max(roots([airinlet ...
        (airinlet+1)*gamma ...
        gamma*gamma*(1-1  ./ MSCtrue)]));
    sigma2REF   = airinlet*sigma2UT;
    % gamma*gamma/(gamma+sigma2UT)/(gamma+sigma2REF)
    
    wREF  = randn(N,1);
    wUT   = randn(N,1);
    
    xREF  = BGfield  + sqrt(sigma2REF)*wREF;
    xUT   = BGfield  + sqrt(sigma2UT)*wUT;
    
    yREF    = filter(bREF,aREF,xREF);
    yUT     = filter(bUT,aUT,xUT);
    % U = 1, REF =2
    signals = [yUT, yREF];
    
    for im=1:LM
        M  = listM(im);
        P=1;
        filtercharacteristics(P).designname         = '';
        filtercharacteristics(P).Norder             = 0;
        filtercharacteristics(P).Wlow_Hz            = 0.0001;
        filtercharacteristics(P).Whigh_Hz           = Fs_Hz;
        filtercharacteristics(P).SCPperiod_sec      = FFT_sec*M;
        filtercharacteristics(P).windowshape        = 'hann';
        filtercharacteristics(P).overlapDFT         = 0.5;
        filtercharacteristics(P).overlapSCP         = 0;
        filtercharacteristics(P).ratioDFT2SCP       = M;
        
        [SUTs, filteredsignals, allfrqsFFT_Hz, alltimes_sec, filterbank] = ...
            fbankanalysis(signals, ...
            filtercharacteristics, ...
            Fs_Hz,...
            0.5);
        
        nbruns = size(SUTs.estimRsup.tabmod,2);
        rmse_integrated = zeros(nbruns,1);
        for irun=1:nbruns
            HestUUUR_irun = SUTs.estimRsup.tabmod(:,irun) .* abs(HREF);
            mseHestUUUR_irun = (HestUUUR_irun(bandave)-abs(HUT(bandave))) .^2;
            mseHestUUUR_irun_rel = mseHestUUUR_irun ./ abs(HUT(bandave));
            rmse_integrated(irun)=sqrt(nanmean(mseHestUUUR_irun));
        end
        RMSEUUUR(iMSC,im) = nanmean(rmse_integrated);
    end
end
%%
figure(4)
leg = cell(LM,1);
for it=1:LM
    leg{it}=sprintf('M = %3i',listM(it));
end
plot(listMSC,RMSEUUUR,'x-')
% hold on
% plot(listMSC,RMSEAB2,'o-')
set(gca,'fontsize',12,'fontname','times')
grid on
hl=legend(leg);
set(hl,'fonts',12,'fontn','times')
xlabel('coherence','fontsize',12,'fontname','times')
ylabel('RMSE on SUT','fontsize',12,'fontname','times')
hold on
plot([0.8 1],0.05*ones(2,1),'--')
hold off
set(gca,'ylim',[0 0.2])
%===========================================
%%
HorizontalSize = 12;
VerticalSize   = 8;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a4');
set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);

set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
set(gcf, 'InvertHardCopy', 'off');

% figure(4); print -depsc -loose ../../textes/6distConjointHMSC/figures/allHest
% print -depsc -loose ../../textes/6distConjointHMSC/figures/allHest
% !epstopdf ../../textes/6distConjointHMSC/figures/allHest.eps
% !rm ../../textes/6distConjointHMSC/figures/allHest.eps

return
%%
% save signals signals Fs_Hz
% return
figure(4)
subplot(411)
semilogx(frqs,10*log10(sRRf))
set(gca,'xlim',[0 6])
set(gca,'ylim',[-20 20])
ylabel('dB')
grid on
%==
subplot(412)
semilogx(frqs,10*log10(sUUf))
set(gca,'xlim',[0 6])
set(gca,'ylim',[-20 20])
ylabel('dB')
grid on
%==
subplot(413)
semilogx(frqs,20*log10(mean(HestUUUR,2)))
hold on
semilogx(frqs,20*log10(abs(HUT)),'r')
hold off
set(gca,'ylim',20*log10(abs(HUT(1)))+[-15,5])
% set(gca,'ylim',[-20 40])
set(gca,'xlim',[0 6])
ylabel('dB')
grid on
%==
subplot(414)
semilogx(frqs,MSC)
set(gca,'xlim',[0 6])
grid on


%
HorizontalSize = 25;
VerticalSize   = 18;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a4');
% set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gca,'fontn','times','fonts',10)

set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');





