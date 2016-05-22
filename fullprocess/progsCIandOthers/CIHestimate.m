%=======================================================================
%   CIHestimate.m
%=======================================================================
% We compar theoretical and simulated distributions
% the thoeretical values are performed in the function
%                 statsRatiosHbis.m (see toolbox)
%=======================================================================
% statistics by drawing vectors under a given spectral matrix GAMMA
% WARNING: that means that the frequency value is fixed and hence can be
% omitted.
%===
% Comparison of Monte-Carlo simulation results and
% theoretical values for the two ratios used to estimate
% HU from HR. The results do depend mainly on the
% true ratio denoted
%
%       HUonHR_k = HU_k / HR_k
%
% which is a complex number.
% To obtain consistent estimates we mix 2M+1 values, either
% by smoothing the periodogram or by averaging several periodograms
% (Welch's approach).
% The statitstics of the modulus of HUonHR_k are given in the report.
% For the phase we used an approximative formula based on delta-method
% and gaussian shape.
%
%%
% Important remark: 
% Below we also compute the conditional bias 
% w.r.t. a given MSC threshold. Therefore when we reduce "threshold2"
% the bias decreases but that does not simulate what we have in 
% real situation. Indeed in this simulaiton the true MSC stays equal 
% to "MSC_true"
% Therefore we are not in the case where we threshold a real MSC 
% which does change.
%=======================================================================

% Reference =1, Under test=2
%       REF UxR
%       UxR UT
%=======================================================================
%
clear
addpath  ../ZZtoolbox/

%====================================================
% R=1, U=2
% we loock for statistics of ratioUUonUR, ratioURonRR,
% then HU = HR x ratioUUonUR OR HU = HR x ratioURonRR
%====================================================

% global factor for the spectral matrices(no effect)
gammafactor            = 5;
absHUonHR              = 1;
phase_HUminusHR_degree = 20;
argHUonHR_rad          = phase_HUminusHR_degree*pi/180;
HRonHU                 = absHUonHR*exp(1j*argHUonHR_rad);
twoMminus1              = 9;
% noise on the sensor under test is currently lower
% than noise on the reference sensor, ie
%      sigmaR = g x sigmaU, with g about 10
% C =1/(1+sigmaU2)(1+sigmaR2)
%
% (1+sigmaU2)(1+g^2 sigmaU2)=1/C
% [1-1/C 1+g^2 g^2]
%%
MSC_true               = 0.96;
Ninlets                = 96;
noiseratioEff          = sqrt(Ninlets);
sigmaU2                = roots([noiseratioEff^2  1+noiseratioEff^2 1-1/MSC_true]);
sigmaU2                = max(sigmaU2);
sigmaU                 = sqrt(sigmaU2);
sigmaR                 = sigmaU*noiseratioEff;
sigmaR2                = sigmaR ^2;
%%
%=========================
%====== spectralmatrix
spectralmatrix         = gammafactor * ...
    [(1+sigmaR2)/absHUonHR  exp(-1j*argHUonHR_rad); ...
    exp(1j*argHUonHR_rad) (1+sigmaU2)*absHUonHR];
%=========================
% for histograms
nbbins = 500;

absHR2 = gammafactor/absHUonHR;
HR     = sqrt(absHR2);
% Monte-Carlo simulations
Lruns  = 1000000;
sqrtG  = sqrtm(spectralmatrix);
SUU_MC = zeros(Lruns,1);
SRR_MC = zeros(Lruns,1);
SUR_MC = zeros(Lruns,1);

for in = 1:twoMminus1
    W  = (randn(Lruns,2)+1j*randn(Lruns,2))/sqrt(2);
    X         = W * sqrtG;
    % REF is on RR=1, and UT is on UU=2
    SRR_MC    = SRR_MC + (X(:,1) .* conj(X(:,1)));
    SUU_MC    = SUU_MC + (X(:,2) .* conj(X(:,2)));
    SUR_MC    = SUR_MC + (X(:,1) .* conj(X(:,2)));
end
SUU_MC           = SUU_MC/twoMminus1;
SRR_MC           = SRR_MC/twoMminus1;
SUR_MC           = SUR_MC/twoMminus1;
absSUR_MC        = abs(SUR_MC);

Rsup_MC          = SUU_MC ./ SUR_MC;
realRsup_MC      = real(Rsup_MC);
imagRsup_MC      = imag(Rsup_MC);
phaseRsup        = -atan2(imagRsup_MC,realRsup_MC);
phaseRsup_degree = phaseRsup*180/pi;

Rinf_MC          = conj(SUR_MC) ./ SRR_MC;
realRinf_MC      = real(Rinf_MC);
imagRinf_MC      = imag(Rinf_MC);
phaseRinf        = -atan2(imagRinf_MC,realRinf_MC);
phaseRinf_degree = phaseRsup*180/pi;

%==================================================
absRsup_MC = abs(Rsup_MC);
absRinf_MC = abs(Rinf_MC);

absRmid_MC = sqrt(SUU_MC ./ SRR_MC);

%============== compare mean/median ================
medianrealRsup         = nanmedian(realRsup_MC);
medianimagRsup         = nanmedian(imagRsup_MC);
absmedianRsup          = sqrt(medianrealRsup .^2 + ...
    medianimagRsup .^2);
phaseRmediansup        = atan2(medianimagRsup,medianrealRsup);
phaseRmediansup_degree = phaseRmediansup*180/pi;

meanrealRsup         = nanmean(realRsup_MC);
meanimagRsup         = nanmean(imagRsup_MC);
absmeanRsup          = sqrt(meanrealRsup .^2 + ...
    meanimagRsup .^2);
phaseRmeansup        = atan2(meanimagRsup,meanrealRsup);
phaseRmeansup_degree = phaseRmeansup*180/pi;
%====================================================

hatsigmau   = sigmaU*(1+randn(Lruns,1)/sqrt(twoMminus1));
hatsigmar   = sigmaR*(1+randn(Lruns,1)/sqrt(twoMminus1));

hatsigmau2  = hatsigmau .^2;
hatsigmar2  = hatsigmar .^2;

AK               = (SUU_MC .* absHR2) ./ (absSUR_MC .^2);
hatgammaSOI      = 0.5*(1+sqrt(1+4*hatsigmau2 .* AK)) ./ AK;
hatabslambda_MC  = absSUR_MC ./ (hatgammaSOI .* absHR2);

%==================================================
noiseratioEff2=noiseratioEff^2;
if noiseratioEff>1
    Hgknown_MC = (((noiseratioEff2-1)/(2*noiseratioEff2)) ...
        .*(absRsup_MC)) .*...
        (1+sqrt(1+4*noiseratioEff2/((noiseratioEff2-1)^2)*...
        absRinf_MC ./ absRsup_MC));
end

MSC_MC = real((absSUR_MC .^2) ./ (SUU_MC .*SRR_MC));

%=========================
% histograms
[hRinf,binRinf] = hist(absRinf_MC,nbbins);
pdfRinf_MC      = hRinf /Lruns/(binRinf(2)-binRinf(1));

[hRsup,binRsup] = hist(absRsup_MC,nbbins);
pdfRsup_MC      = hRsup /Lruns/(binRsup(2)-binRsup(1));

[hRmid,binRmid] = hist(absRmid_MC,nbbins);
pdfRmid_MC      = hRmid /Lruns/(binRmid(2)-binRmid(1));

[harg,binarg]  = hist(phaseRsup_degree,nbbins);
pdfarg_MC      = harg /Lruns/(binarg(2)-binarg(1));

if noiseratioEff>1
    [hgknown,bingknown] = hist(Hgknown_MC,nbbins);
    pdfgknown_MC        = hgknown /Lruns/(bingknown(2)-bingknown(1));
end

mean2expH        = -(absRinf_MC - absRsup_MC);
[hmean,binmean]  = hist(mean2expH,nbbins);
pdfmean_MC       = hmean /Lruns/(binmean(2)-binmean(1));

[hmym,binmym]    = hist(hatabslambda_MC,100);
pdfmym_MC        = hmym /Lruns/(binmym(2)-binmym(1));

%============ variable ranges
allT.TUUonUR     = binRsup;
allT.TURonRR     = binRinf;
allT.MSC         = linspace(0.6,1,100);
allT.phase       = binarg;
[hatpdfGsup, hatpdfGinf, hatMSC, hatPhase] = ...
    theoreticalStats(allT,spectralmatrix,twoMminus1,0.05);
[hatpdfGsupP, hatpdfGinfP, hatMSCP, hatPhaseP] = ...
    theoreticalStatsV2(allT,spectralmatrix,twoMminus1,0.05);

%%
stdapMLEdeltamethod                 = ...
    sqrt(absHUonHR * absHUonHR * (1-MSC_true)/sqrt(MSC_true)/(2*twoMminus1));
hatpdfappMLE                        = ...
    (1/sqrt(2*pi)/stdapMLEdeltamethod)* ...
    exp(-(binmym-1) .^2/(2 * stdapMLEdeltamethod^2));
%=========== conditional expectation
% E(Rsup|C\in(c1,c2))
LLlistind = 200;
listind   = linspace(0.1,1,LLlistind);
Econd     = zeros(LLlistind-1,1);
for ib=2:LLlistind
    condindex   = and(MSC_MC>listind(ib-1), MSC_MC<listind(ib));
    Econd(ib-1) = nanmean(absRsup_MC(condindex));
end
%%
if 1
    figure(1)
    clf
    subplot(2,1,1)
    bar(binRsup,pdfRsup_MC)
    hold on
    plot(binRsup,hatpdfGsup.pdf,'r','linew',2)
    hold off
   
    figure(2)
    clf
    subplot(2,1,1)
    bar(binRsup,pdfRsup_MC)
    hold on
    plot(binRsup,hatpdfGsupP.pdf,'r','linew',2)
    hold off

    %== title
    title(sprintf('sensor gain ratio = %5.1f, true MSC = % 4.2f,\nnoise ratio = %2i, M = %i',...
        absHUonHR, MSC_true, Ninlets, (twoMminus1+1)/2),'fontsize',12)
    
   %===========================================
    
       figure(1)
    subplot(2,1,2)
    bar(binRinf,pdfRinf_MC)
    hold on
    plot(binRinf,hatpdfGinf.pdf,'r','linew',2)
    hold off
    
    figure(2)
    subplot(2,1,2)
    bar(binRinf,pdfRinf_MC)
    hold on
    plot(binRinf,hatpdfGinfP.pdf,'r','linew',2)
    hold off
    
    %     %===========================================
    %
    %     subplot(3,1,3)
    %     %     bar(bingknown,pdfgknown_MC)
    %     %     set(gca,'xlim',[0.2 1.8])
    %     bar(binarg,pdfarg_MC)
    %     hold on
    %     plot(binarg,hatPhase,'r','linew',2)
    %     hold off
    %     grid on
    %     % text('string',txt,'interpreter','latex','pos',[0.75,3],'fontsize',16)
    %     set(gca,'fontname','times','fontsize',12)
    %===========================================
    txtsupinf    = cell(2,1);
%     textEsup = sprintf('mean-sup = %4.5f', nanmean(absRsup_MC));
%     textEinf = sprintf('mean-inf = %4.5f', nanmean(absRinf_MC));
    textEsup = sprintf('mean-sup = %4.3f', hatpdfGsup.mean);
    textEinf = sprintf('mean-inf = %4.3f', hatpdfGinf.mean);

    txtsupinf{1} = '$\widehat{R}_{\sup}$';
    txtsupinf{2} = '$\widehat{R}_{\inf}$';
    
    textR{1} = sprintf('%s\n%s',txtsupinf{1},textEsup);
    textR{2} = sprintf('%s\n%s',txtsupinf{2},textEinf);

for ii=1:2
    subplot(2,1,ii)
    set(gca,'fontname','times','fontsize',12)
    text('string',textR{ii},...
        'interpreter','latex','pos',[1.08 4.8],'fontsize',12)
    set(gca,'xlim',nanmean(absRsup_MC) + 0.3*[-1 1])
    set(gca,'ylim',[0 10])
    grid on
end

    HorizontalSize = 12;
    VerticalSize   = 8;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a4');
%     set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    
    set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
    set(gcf, 'InvertHardCopy', 'off');
end

% figure(1); print -depsc -loose ../../textes/6distConjointHMSC/figures/theoreticaldistribratios.eps

%==========================================
if 0
    threshold1             = MSC_true;
    threshold2             = 0.98;
    figure(2)
    clf
    subplot(121)
    plot(MSC_MC(1:Lruns),(absRsup_MC(1:Lruns)),'.')
    set(gca,'ylim',[0 2.5])
    set(gca,'xlim',[0 1])
    grid on
    hold on
    plot([0,1],absHUonHR*ones(2,1),'k','linew',2)
    hold off
    set(gca,'fontname','times','fontsize',18)
    xlabel('MSC','fontsize',18)
    ylabel('Rsup','fontsize',18)
    
    % subplot(223)
    % plot(listind(2:LLlistind)-0.5* (listind(2)-listind(1)),...
    %     20*log10(Econd))
    % grid on
    % set(gca,'fontn','times','fonts',10)
    % ylabel('Conditional expectation - dB')
    
    text00 = sprintf('Simulation based on spectral matrix');
    text0 = sprintf('True ratio HUT/HREF = %i,\nTrue MSC = % 4.2f\nNoise power ratio = %2i (%4.1f dB)\nSmoothing window number = %i\nLruns = %i',  ...
        absHUonHR,MSC_true, Ninlets, 20*log10(noiseratioEff), twoMminus1, Lruns);
    text1 = sprintf('E(Rsup|MSC>%4.2f) - 1 = %4.5f', threshold1, ...
        nanmean(absRsup_MC(MSC_MC>threshold1))-1);
    text2 = sprintf('E(Rsup|MSC>%4.2f) - 1 = %4.5f', threshold2, ...
        nanmean(absRsup_MC(MSC_MC>threshold2))-1);%, sum(MSC_MC>threshold2));
%     text3 = sprintf('E(Rsup) = %4.5f', nanmean(absRsup_MC));
    text4 = sprintf('median(Rsup) = %4.5f', nanmedian(absRsup_MC));
%     text5 = sprintf('E(Rinf) = %4.5f', nanmean(absRinf_MC));
    txts = sprintf('%s\n%s\n%s\n%s\n%s\n%s',text00,text0,text1,text2,textEsup,textEinf);
    
    subplot(122)
    set(gca,'ylim',[0 0.2])
    set(gca,'xtick',[],'ytick',[],'box','off')
    set(gca,'color', [1,1,0.92])
    set(gca,'xcolor', [1,1,0.92])
    set(gca,'ycolor', [1,1,0.92])
    text(0,0.1,txts,'fontname','times','fontsize',18)
    hpos = get(gca,'position');
    set(gca,'position',[0.52 hpos(2:4)])
    
    set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
    set(gcf, 'InvertHardCopy', 'off');

end
%%
printdirectory  = ' ../../textes/6distConjointHMSC/slidesITW2015/';
fileprintepscmd = sprintf('print -depsc -loose %stheoreticaldistribratios.eps',printdirectory);
fileeps2pdfcmd  = sprintf('/Library/TeX/texbin/epstopdf %stheoreticaldistribratios.eps',printdirectory);
filermcmd       = sprintf('!rm %stheoreticaldistribratios.eps',printdirectory);
saveflag=0;
%
if saveflag
    eval(fileprintepscmd)
%     eval(fileeps2pdfcmd)
    eval(filermcmd)
end

