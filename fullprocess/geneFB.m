clear

allcolors = ['g.';'y.';'m.';'r.';'k.';'b.';'rx';'yx';'mx';'rx';'kx';'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];

M                               = 5;

filtercharact = struct;
filtercharact(1).SCPperiod_sec  = 1000 ;
filtercharact(2).SCPperiod_sec  = 500 ;
filtercharact(3).SCPperiod_sec  = 250 ;
filtercharact(4).SCPperiod_sec  = 100 ;
filtercharact(5).SCPperiod_sec  = 50 ;
filtercharact(6).SCPperiod_sec  = 25;

filtercharact(1).overlapDFT = 0.5;
filtercharact(2).overlapDFT = 0.5;
filtercharact(3).overlapDFT = 0.5;
filtercharact(4).overlapDFT = 0.5;
filtercharact(5).overlapDFT = 0.5;
filtercharact(6).overlapDFT = 0.5;

filtercharact(1).overlapSCP = 0;
filtercharact(2).overlapSCP = 0;
filtercharact(3).overlapSCP = 0;
filtercharact(4).overlapSCP = 0;
filtercharact(5).overlapSCP = 0;
filtercharact(6).overlapSCP = 0;

filtercharact(1).ratioDFT2SCP = M;
filtercharact(2).ratioDFT2SCP = M;
filtercharact(3).ratioDFT2SCP = M;
filtercharact(4).ratioDFT2SCP = M;
filtercharact(5).ratioDFT2SCP = M;
filtercharact(6).ratioDFT2SCP = M;

filtercharact(1).Norder  = 2 ;
filtercharact(2).Norder  = 3 ;
filtercharact(3).Norder  = 4 ;
filtercharact(4).Norder  = 4 ;
filtercharact(5).Norder  = 4 ;
filtercharact(6).Norder  = 6 ;

filtercharact(1).designname     = 'butter';
filtercharact(1).windowshape    = 'hann';
filtercharact(1).Wlow_Hz        = 0.005;
filtercharact(1).overlapDFT     = 0.5;
TFFT_Pfilter                    = filtercharact(1).SCPperiod_sec/M;

for Pfilter=2:6
    filtercharact(Pfilter).designname     = 'butter';
    filtercharact(Pfilter).windowshape    = 'hann';
    filtercharact(Pfilter).overlapDFT     = 0.5;
%     filtercharact(Pfilter).SCPperiod_sec  = filtercharact(Pfilter-1).SCPperiod_sec*0.5;
    TFFT_Pfilter                          = filtercharact(Pfilter).SCPperiod_sec/M;
    filtercharact(Pfilter).Wlow_Hz        = (1/TFFT_Pfilter)/0.2;
     filtercharact(Pfilter-1).Whigh_Hz     = 1.1*filtercharact(Pfilter).Wlow_Hz;
    filtercharact(Pfilter).Wlow_Hz        = (1/TFFT_Pfilter)/0.08;
     filtercharact(Pfilter-1).Whigh_Hz     = filtercharact(Pfilter).Wlow_Hz;
    filtercharact(Pfilter).Wlow_Hz        = (1/TFFT_Pfilter)/0.2;
     filtercharact(Pfilter-1).Whigh_Hz     = 1.1*filtercharact(Pfilter).Wlow_Hz;
end
filtercharact(Pfilter).Whigh_Hz           = 8;

scal1 = 1 ./ ([filtercharact(:).SCPperiod_sec] .* ...
    ([filtercharact(:).Whigh_Hz]-[filtercharact(:).Wlow_Hz]));

scal2 = 2 ./ ([filtercharact(:).SCPperiod_sec] .* ...
    ([filtercharact(:).Whigh_Hz]+[filtercharact(:).Wlow_Hz]));

Fs_Hz                 = 20;
P                     = length(filtercharact);
Lfft                  = 16*2048;
frqs_Hz               = (0:Lfft-1)*Fs_Hz/Lfft;
%========== plots of the response
figure(1)
clf
%========== end

for ip=1:P
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
for ip=1:P
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
for ip = 1:P
    tabFB = [tabFB ...
        sprintf('$[%4.3f-%4.2f]$&$%i$\n%s %s',...
        filtercharact(ip).Wlow_Hz, ...
        filtercharact(ip).Whigh_Hz,filtercharact(ip).SCPperiod_sec ,'\\ \hline')];
end
tabFB
%====
save filtercharacteristics2 filtercharact

filtercharact
