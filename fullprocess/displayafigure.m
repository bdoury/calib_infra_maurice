%%
figdisplay = 3;
figure(figdisplay)
%==============================================
allcolors = ['b.';'r.';'m.';'c.';'g.';'k.';'rx';'yx';'mx';'rx';'kx';...
    'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];
if and(P==1,isnan(filterbank{1}.num))
    ip = 1;
    fq_ip = SUTs(ip).frqsFFT_Hz;
    Rsup_ip = 20*log10(SUTs(ip).estimRsup.modcst);
    Rinf_ip = 20*log10(SUTs(ip).estimRinf.modcst);
    racN_ip = sqrt(SUTs(ip).Nsupthresholdintheband);
    stdRsup_ip = SUTs(ip).estimRsup.stdmodcst ./racN_ip ;
    semilogx(fq_ip,Rsup_ip,...
        'o','color',allcolors(ip,1),'markerfac',allcolors(ip,1))
    hold on
    %             semilogx(fq_ip,Rinf_ip,...
    %                 'o','color',allcolors(ip,1),'markerfac','k')
    semilogx(fq_ip, Rsup_ip+2*stdRsup_ip,'o','color',allcolors(ip,1))
    semilogx(fq_ip, Rsup_ip-2*stdRsup_ip,'o','color',allcolors(ip,1))
    semilogx(ones(2,1)*fq_ip, ...
        [Rsup_ip-2*stdRsup_ip,Rsup_ip+2*stdRsup_ip]','color',allcolors(ip,1),'linew',1.5)
else
    idipinf = zeros(P,1);
    idipsup = zeros(P,1);
    for ip=1:P
        idipinf(ip) = SUTs(ip).indexinsidefreqband(1);
        idipsup(ip) = SUTs(ip).indexinsidefreqband(2);
    end
    for ip=1:P
        fq_ip = SUTs(ip).frqsFFT_Hz(idipinf(ip):idipsup(ip));
        Rsuptab_ip = 20*log10(SUTs(ip).estimRsup.tabmodcst(idipinf(ip):idipsup(ip),:));
        semilogx(fq_ip, Rsuptab_ip,'.','color',[1 1 1]*0.9)
    end
    for ip=1:P
        fq_ip = SUTs(ip).frqsFFT_Hz(idipinf(ip):idipsup(ip));
        Rsup_ip = 20*log10(SUTs(ip).estimRsup.modcst(idipinf(ip):idipsup(ip)));
        Rinf_ip = 20*log10(SUTs(ip).estimRinf.modcst(idipinf(ip):idipsup(ip)));
        racN_ip = sqrt(SUTs(ip).Nsupthresholdintheband);
        stdRsup_ip = SUTs(ip).estimRsup.stdmodcst(idipinf(ip):idipsup(ip)) ./racN_ip ;
        semilogx(fq_ip,Rsup_ip,...
            'o','color',allcolors(ip,1),'markerfac',allcolors(ip,1))
        hold on
        %             semilogx(fq_ip,Rinf_ip,...
        %                 'o','color',allcolors(ip,1),'markerfac','k')
        semilogx(fq_ip, Rsup_ip+2*stdRsup_ip,'o','color',allcolors(ip,1))
        semilogx(fq_ip, Rsup_ip-2*stdRsup_ip,'o','color',allcolors(ip,1))
        semilogx(ones(2,1)*fq_ip, ...
            [Rsup_ip-2*stdRsup_ip,Rsup_ip+2*stdRsup_ip]','color',allcolors(ip,1),'linew',1.5)
    end
end
hold off
grid on
set(gca,'xlim',[1e-3 6])
set(gca,'ylim',[-2 2])
xlabel('frequency - Hz')
ylabel('gain ratio - dB')
set(gca,'fontname','times','fontsize',16)
%====================
figure(figdisplay)
% subplot(121)
hold off
title('Infrasonic signals - IS26','fontsize',24)
% subplot(122)
% hold off
HorizontalSize = 14;
VerticalSize   = 7;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
%         set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');
