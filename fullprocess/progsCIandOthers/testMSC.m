clear all
% warning: for CI we have to take a/2 and 1-a/2
% for monolateral test we have to take 1-a
%
addpath  /Users/maurice/etudes/ctbto/allJOBs2015/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/

listC0      = 0.0:0.01:1;
listNd      = 5:1:8;
LC          = length(listC0);
LN          = length(listNd);
a           = 0.1;
CI          = zeros(LC,3,LN);
for iN = 1:LN
    Nd = listNd(iN)-1;
    for ic = 1: LC
        C0=listC0(ic);
        CI(ic,1,iN)=invcumulFunctionMSC(a/2,C0,Nd);
        CI(ic,2,iN)=invcumulFunctionMSC(1-a/2,C0,Nd);
        CI(ic,3,iN)=invcumulFunctionMSC(1-a,C0,Nd);
    end
end
%%
figure(1)
subplot(121)
plot(listC0,squeeze(CI(:,3,:)))
hold on
plot([0,1],[0,1],'--')
hold off
grid on
ylabel('\eta')
xlabel('C')
set(gca,'xlim',[0  1])
set(gca,'ylim',[0  1])

subplot(122)
plot(listC0,squeeze(CI(:,3,:)))
set(gca,'xlim',[0.85 0.95])
set(gca,'ylim',[0.95 1])
hold on
hplot = plot(0.9,squeeze(CI(listC0==0.9,3,:)),'o','markersize',8);
for iN=1:LN
    set(hplot(iN),'markerf',get(hplot(iN),'color'))
end
hold off
grid on
ylabel('\eta')
xlabel('C')
legend(sprintf('N_d = %i',listNd(1)), ...
    sprintf('N_d = %i',listNd(2)), ...
    sprintf('N_d = %i',listNd(3)), ...
    sprintf('N_d = %i',listNd(4)))

HorizontalSize = 18;
VerticalSize   = 12;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a4');
% set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gca,'fontname','times','fontsize',10)

set(gcf,'color', [1,1,0.92]);%0.7*ones(3,1))
set(gcf, 'InvertHardCopy', 'off');

%  print -depsc -loose  ../texte2/thresholdforC95.eps

%===
return
listthreshold = (0:0.05:1)';
LTh           = length(listthreshold);
ROC           = zeros(LTh, LN,LC,2);

for iN = 1:LN
    Nd = listNd(iN);
    for ic = 1: LC
        C0=listC0(ic);
        ROC(:,iN,ic,1) = cumulFunctionMSC(listthreshold,0,Nd);
        ROC(:,iN,ic,2) = cumulFunctionMSC(listthreshold,C0,Nd);
    end
end
