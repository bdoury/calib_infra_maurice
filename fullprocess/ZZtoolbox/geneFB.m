%============== geneFB.m ==========================
% this program built the file .mat which consists  
% of thecharacteristic of the time/frequency 
% pavement.
% the array of the following structures is saved 
% in a file .mat
%==================================================
%     SCPperiod_sec
%     overlapDFT: typically 0.5
%     overlapSCP: typically 0
%     ratioDFT2SCP: 
%     Norder
%     designname
%     windowshape
%     Wlow_Hz
%     Whigh_Hz
%==================================================
% the length in second of the FFT is :
%           SCPperiod_sec/ratioDFT2SCP
% If overlapDFT=0.5, and if ratioDFT2SCP=5,
% the number of FFTs to perform the SCP is 9
% If overlapSCP=0, the next SCP is performed by shifting
% the FFT block with 9
%
%===================================================

clear

commandsave = 'save filtercharacteristics/filtercharacteristics filtercharact';

allcolors = ['g.';'y.';'m.';'r.';'k.';'b.'; ...
    'rx';'yx';'mx';'rx';'kx';'c.';'k.';'r.';'c.';'m.';'g.';...
    'b.';'k.';'r.';'c.';'m.';'g.';'k.'];

M                               = 5;
Pfilter                         = 6;
listofperiods                   = [1000,800,350,150,80,30];
listOrder                       = [2 2 2 2 3 4];
filtercharact = struct;

%==== Remark =====================================
% due to some constraints in programming
% the value filtercharact(iP).overlapSCP must be
% 0, 1/M, 2/M, 3/M, ..., M-1/M which corresponds to
% the frontier of the FFT block.
% In pratice we advise to take 0
    
for iP = 1:Pfilter
    filtercharact(iP).SCPperiod_sec  = listofperiods(iP) ;
    filtercharact(iP).overlapDFT = 0.5;
    filtercharact(iP).overlapSCP = 0;
    filtercharact(iP).ratioDFT2SCP = M;
    filtercharact(iP).Norder  = listOrder(iP) ;
    filtercharact(iP).designname     = 'butter';
    filtercharact(iP).windowshape    = 'hann';
    filtercharact(iP).overlapDFT     = 0.5;
end
filtercharact(1).Wlow_Hz        = 0.005;
filtercharact(1).overlapDFT     = 0.5;
TFFT_Pfilter                    = filtercharact(1).SCPperiod_sec/M;
for iP=2:Pfilter
    %     filtercharact(Pfilter).SCPperiod_sec  = filtercharact(Pfilter-1).SCPperiod_sec*0.5;
    TFFT_Pfilter                          = filtercharact(iP).SCPperiod_sec/M;
    filtercharact(iP).Wlow_Hz        = (1/TFFT_Pfilter)/0.2;
    filtercharact(iP-1).Whigh_Hz     = 1.1*filtercharact(iP).Wlow_Hz;
    filtercharact(iP).Wlow_Hz        = (1/TFFT_Pfilter)/0.08;
    filtercharact(iP-1).Whigh_Hz     = filtercharact(iP).Wlow_Hz;
    filtercharact(iP).Wlow_Hz        = (1/TFFT_Pfilter)/0.2;
    filtercharact(iP-1).Whigh_Hz     = 1.1*filtercharact(iP).Wlow_Hz;
end
filtercharact(iP).Whigh_Hz           = 8;

scal1 = 1 ./ ([filtercharact(:).SCPperiod_sec] .* ...
    ([filtercharact(:).Whigh_Hz]-[filtercharact(:).Wlow_Hz]));

scal2 = 2 ./ ([filtercharact(:).SCPperiod_sec] .* ...
    ([filtercharact(:).Whigh_Hz]+[filtercharact(:).Wlow_Hz]));

Fs_Hz                 = 20;
Lfft                  = 16*2048;
frqs_Hz               = (0:Lfft-1)*Fs_Hz/Lfft;
%========== plots of the response
figure(1)
clf
%========== end

for ip=1:Pfilter
    subplot(121)
    [filnum,filden] = butter(filtercharact(ip).Norder,...
        2*[filtercharact(ip).Wlow_Hz ...
        filtercharact(ip).Whigh_Hz]/Fs_Hz);
    Hf = abs(fft(filnum,Lfft) ./ fft(filden,Lfft));
    semilogx(frqs_Hz, 10*log10(Hf),allcolors(ip),'linew',2)
    hold on
    set(gca,'xlim',[1e-3 10])
    set(gca,'ylim',[-20 2])
    set(gca,'color',[1 1 1]*0.9)
    grid on
end
hold off
set(gca,'fontname','times','fontsize',12);
xlabel('frequency - Hz')
ylabel('gain - dB')
loc = get(gca,'position');
% bar graphics, length proportional to the duration
loc(1)  = loc(1)+0.4;
loc(2)  = 0.92;
highbar = 0.12;
for ip=1:Pfilter
    switch ip
        case 1
            lenbar  = (filtercharact(ip).SCPperiod_sec)/2200;
            locx = loc(1)+loc(3)*0.12-0.14;
        case 2
            lenbar  = (filtercharact(ip).SCPperiod_sec)/2200;
            locx = loc(1)+loc(3)*0.12;
        case 3
            lenbar  = (filtercharact(ip).SCPperiod_sec)/1800;
            locx = loc(1)+loc(3)*0.12+0.6;
        case 4
            lenbar  = (filtercharact(ip).SCPperiod_sec)/1700;
            locx = loc(1)+loc(3)*0.12+0.7;
        case 5
            lenbar  = (filtercharact(ip).SCPperiod_sec)/1700;
            locx = loc(1)+loc(3)*0.12+1;
        case 6
            lenbar  = (filtercharact(ip).SCPperiod_sec)/1700;
            locx = loc(1)+loc(3)*0.12+1;
    end
    midy    = (filtercharact(ip).Whigh_Hz+filtercharact(ip).Wlow_Hz)/12;
    loc(2)  = loc(2)-highbar;
    subplot('position',[loc(1) loc(2) lenbar highbar])
    set(gca,'xtick',[],'ytick',[],'box','on')
    tt=sprintf('%i seconds - [%4.2f %4.2f] Hz',filtercharact(ip).SCPperiod_sec,...
        filtercharact(ip).Wlow_Hz,...
        filtercharact(ip).Whigh_Hz);
    text(locx,midy*0.03+0.355,tt,'fontname','times','fontsize',14)
    set(gca,'color',allcolors(ip))
end
set(gca,'fontname','times','fontsize',12);
% xlabel('time of stationarity - frequency bandwidth')
HorizontalSize = 26;
VerticalSize   = 8;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
% set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);

set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');

%
%         print -depsc -loose ../../textes/6distConjointHMSC/filterbank.eps
%     !epstopdf ../../textes/6distConjointHMSC/filterbank.eps
%     !rm ../../textes/6distConjointHMSC/filterbank.eps
%      print -depsc -loose ../../slides/7sessionmeeting/filterbank.eps
%         print -depsc -loose ../../slides/6slidessummary01072915/.eps
%     !epstopdf ../../slides/6slidessummary01072915/filterbank.eps
%     !rm ../../slides/6slidessummary01072915/filterbank.eps
%
%     printdirectory  = ' ../slidesITW2015/';
%     fileprintepscmd = sprintf('print -depsc -loose %sfilterbank.eps',printdirectory);
%     fileeps2pdfcmd  = sprintf('!epstopdf %sfilterbank.eps',printdirectory);
%     filermcmd       = sprintf('!rm %sfilterbank.eps',printdirectory);
%     %
%       eval(fileprintepscmd)
%         eval(fileeps2pdfcmd)
%         eval(filermcmd)


%==== parameters for TeX
tabFB = [];
for ip = 1:Pfilter
    tabFB = [tabFB ...
        sprintf('$[%4.3f-%4.2f]$&$%i$\n%s %s',...
        filtercharact(ip).Wlow_Hz, ...
        filtercharact(ip).Whigh_Hz,filtercharact(ip).SCPperiod_sec ,'\\ \hline')];
end
tabFB
%====
eval(commandsave)

filtercharact
