%===============================================================
% This program draws the frequency responses of
% the filter bank described in the file 
%      'filtercharactfilename'
% Also are represented the times in second
% associated to each filters.
% 
% Specifically the filter element writes:
% %DFT means discrete fourier transform
% %SCP means spectral components
% P=1;
% filtercharact(P).designname      = 'butter';
% filtercharact(P).Norder          = 2;
% filtercharact(P).Wlow_Hz         = 0.02;
% filtercharact(P).Whigh_Hz        = 0.08;
% filtercharact(P).SCPperiod_sec   = 500;
% filtercharact(P).windowshape     = 'hann';
% filtercharact(P).overlapDFT      = 0.5;
% filtercharact(P).overlapSCP      = 0;
% filtercharact(P).ratioDFT2SCP    = 5;
% %------
% P=P+1; etc
% 
%===============================================================
clear all
figure(1)
% close all
filtercharactfilename = 'filtercharacteristics';
allcolors = ['g.';'y.';'m.';'r.';'g.';'c.';'rx';'yx';'mx';'rx';'kx';'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];

%===============================================================
% sampling frequency
Fs_Hz                 = 20;
cmdloadcharact        = sprintf('run(''%s'')',filtercharactfilename);
eval(cmdloadcharact);
P                     = length(filtercharact);
Lfft                  = 16*2048;
frqs_Hz               = (0:Lfft-1)*Fs_Hz/Lfft;
% plots of the response
clf
for ip=1:P
    subplot(121)
    [filnum,filden] = butter(filtercharact(ip).Norder,...
        2*[filtercharact(ip).Wlow_Hz ...
        filtercharact(ip).Whigh_Hz]/Fs_Hz);
    Hf = abs(fft(filnum,Lfft) ./ fft(filden,Lfft));
    semilogx(frqs_Hz, 10*log10(Hf),allcolors(ip),'linew',2)
    hold on
    set(gca,'xlim',[1e-3 10])
    set(gca,'ylim',[-30 2])
    set(gca,'color',[1 1 1]*0.7)
    grid off
end
hold off
set(gca,'fontname','times','fontsize',12);
xlabel('frequency - Hz')
ylabel('gain - dB')
loc = get(gca,'position');
% bar graphics, length proportional to the duration
loc(1) = loc(1)+0.4;
loc(2) = 0.92;
highbar = 0.085;
for ip=1:P
    if (ip==1)
    lenbar  = (filtercharact(ip).SCPperiod_sec)/2600;
    elseif ip==2
            lenbar  = (filtercharact(ip).SCPperiod_sec)/2100;
    else
    lenbar  = (filtercharact(ip).SCPperiod_sec)/1100;
    end
    midy    = (filtercharact(ip).Whigh_Hz+filtercharact(ip).Wlow_Hz)/12;
    loc(2)  = loc(2)-highbar;
    subplot('position',[loc(1) loc(2) lenbar highbar])
    set(gca,'xtick',[],'ytick',[],'box','on')
    tt=sprintf('%i seconds - [%4.2f %4.2f] Hz',filtercharact(ip).SCPperiod_sec,...
        filtercharact(ip).Wlow_Hz,...
        filtercharact(ip).Whigh_Hz);
    text(loc(1)+loc(3)*0.02,midy*0.2+0.5,tt,'fontname','times','fontsize',10)
    set(gca,'color',allcolors(ip))
end
set(gca,'fontname','times','fontsize',12);
% xlabel('time of stationarity - frequency bandwidth')
HorizontalSize = 26;
VerticalSize   = 8;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);

set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');

% 
%         print -depsc -loose ../../textes/6distConjointHMSC/filterbank.eps
%     !epstopdf ../../textes/6distConjointHMSC/filterbank.eps    
%     !rm ../../textes/6distConjointHMSC/filterbank.eps
%      print -depsc -loose ../../slides/7sessionmeeting/filterbank.eps
%         print -depsc -loose ../../slides/6slidessummary01072915/filterbank.eps
%     !epstopdf ../../slides/6slidessummary01072915/filterbank.eps   
%     !rm ../../slides/6slidessummary01072915/filterbank.eps
% 
%==== parameters for TeX
tabFB = [];
for ip = 1:P
    tabFB = [tabFB ...
        sprintf('$[%4.2f-%4.2f]$&$%i$\n%s %s',...
        filtercharact(ip).Wlow_Hz, ...
        filtercharact(ip).Whigh_Hz,filtercharact(ip).SCPperiod_sec ,'\\ \hline')];
end
tabFB
%====


