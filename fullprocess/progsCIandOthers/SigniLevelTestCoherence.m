%==================================================
% SigniLevelTestCoherence.m
%==================================================
% Performances of test on the hypothesis
% H0 = { MCS > cintest} with M d.o.f
% The significance level is denoted signifiancelevel
% Used function: cumulFunctionMSC.m
%==================================================
clear
addpath  /Users/maurice/etudes/ctbto/allJOBs2015/myjob/1TaskOnSensors/textes/6distConjointHMSC/fullprocess/ZZtoolbox/
signifiancelevel = 0.9;
cintest = 0.97; NE=100;
E = linspace((0.9),(0.999),NE)';
allM = 5:2:10;
LM = length(allM);
P = zeros(NE,LM);
thres = zeros(LM,1);
leg=[];
for id=1:LM
    N = allM(id);
    P(:,id) = cumulFunctionMSC(E,cintest,N);
    thres(id) = invcumulFunctionMSC(signifiancelevel,cintest,N);
    leg = [leg;{sprintf('M = %2i, th = %4.3f',N,thres(id))}];
end
plot(E,P,'.-')
hold on
plot([0, 1],[signifiancelevel * [1 1]])
ht = plot(thres,signifiancelevel*ones(LM,1),'o');
set(ht,'markerfacec','b')
hold off
grid
set(gca,'ylim',[0.8 1])
set(gca,'xlim',[0.97 1])

ht=legend(leg);
set(ht,'fontsize',12)
ht=title(sprintf('H_0 = %s MSC>%4.2f %s at level %4.2f','\{',cintest,'\}',signifiancelevel));
set(ht,'fontname','times','fontsize',12)

HorizontalSize = 12;
VerticalSize   = 8;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a4');
set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gca,'fontname','times','fontsize',12)

set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');
%

% figure(4); print -depsc -loose  ../../textes/6distConjointHMSC/figures/MSCtestthreshold

% figure(4)
% print -dpdf -loose ../../textes/6distConjointHMSC/figures/MSCtestthreshold
% print -depsc -loose ../../textes/6distConjointHMSC/figures/allHest
% !epstopdf ../../textes/6distConjointHMSC/figures/allHest.eps
% !rm ../../textes/6distConjointHMSC/figures/allHest.eps
