clear
dirmat = dir('./st*.mat');
nbfiles = length(dirmat);

for ii=1:1%nbfiles
    cmdload = sprintf('load %s',dirmat(ii).name);
    eval(cmdload)

    figure(1)
%     subplot(1,nbfiles,ii)
    leg = cell(L2Mp1,1);
    for it=1:L2Mp1
        leg{it}=sprintf('N = %3i',listtwoMplus1(it));
    end
    plot(listMSC_theo,stdonHUonHRex1,'o-')
    set(gca,'fontsize',16,'fontname','times')
    grid on
    title(sprintf('Gain ratio = %4.2f',absHUonHR))
    hl=legend(leg);
    set(hl,'fontsize',16,'fontname','times')
    xlabel('coherence','fontsize',16,'fontname','times')
    if ii==1
        ylabel('RMSE on SUT','fontsize',16,'fontname','times')
    end
    hold on
    plot([0.8 1],0.05*ones(2,1),'--')
    hold off
    set(gca,'ylim',[0 0.2])
    %===========================================
    HorizontalSize = 21;
    VerticalSize   = 12;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a4');
    set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    
    set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
    set(gcf, 'InvertHardCopy', 'off');
    
%     cmdprint = sprintf('print -dpdf %s.pdf',dirmat(ii).name);
%     eval(cmdprint)
end

cmdprint = sprintf('print -dpdf allHest.pdf');
eval(cmdprint)

