clear

dataresults='AAresultswithFBXsepar/';
indexofSTA = 1;
filesindir = dir(sprintf('%ss%i/*.mat',dataresults,indexofSTA));
nbofcouplesdays = length(filesindir);
allsRatio = zeros(10000,nbofcouplesdays);
for inday=1:nbofcouplesdays
    filename = filesindir(inday).name;
    comload = sprintf('load %ss%i/%s',dataresults,indexofSTA,filename);
    eval(comload);
    if inday==1
        frqs = allfrqsPfilters;
        allsRatio = allsRatio(1:length(frqs),:);
    end
    allsRatio(:,inday) = allRatioPfilters(:,inday);
end
allsRatio             = allsRatio(1:length(frqs),:);
[frqsU, indunique]    = unique(frqs);
allRatioU             = allsRatio(indunique,:);
frqsUZ               = frqsU(not(frqsU==0));
allRatioUZ            = allRatioU(not(frqsU==0),:);
meanAllratioUZ        = nanmean(allRatioUZ,2);
absmeanAllratioUZ     = abs(meanAllratioUZ);
phasemeanAllratioUZ  = angle(meanAllratioUZ);

            segmentsnumber = 2;
            polydegrees    = (0:1:8);
            logtrain.flag  = 1;
            logtrain.N     = 70;
            logfit.flag    = 1;
            logfit.N       = 60;

            [FreqFitabs_Hz,absRfit] = smoothpolyLL(frqsUZ, ...
                absmeanAllratioUZ,...
                segmentsnumber,polydegrees,logtrain,logfit);

        

% figure(1)
% 
% Npermut=10;
% P = randperm(nbofcouplesdays,Npermut);
% meanP = abs(nanmean(allRatioU(:,P),2));
% Q = setdiff((1:nbofcouplesdays),P);
%     semilogx(frqsU,20*log10(meanP),'r')
% hold on
% semilogx(frqsU,20*log10(meanP*1.05),'r')
% semilogx(frqsU,20*log10(meanP*0.95),'r')
% 
% for inday=1:nbofcouplesdays-Npermut
%     semilogx(frqsU,20*log10(abs(allRatioU(:,Q(inday)))),'.','color',0.5*[1 1 1])
% end
% hold off
% grid on

figure(2)
PP = randperm(nbofcouplesdays);
allRatioUZP = allRatioUZ(:,PP);

Nfollowing = 8;
% semilogx(frqsUZ,20*log10(absmeanAllratioUZ),'r')
% hold on
semilogx(frqsUZ,20*log10(absmeanAllratioUZ*1.05),'r')
hold on
semilogx(frqsUZ,20*log10(absmeanAllratioUZ*0.95),'r')

for ii=1:fix(nbofcouplesdays/Nfollowing)
    id1=(ii-1)*Nfollowing+1;
    id2=min([id1+Nfollowing,nbofcouplesdays]) ;
    semilogx(frqsUZ,20*log10(abs(nanmean(allRatioUZP(:,id1:id2),2))),'.','color',0.7*[1 1 1])
    hold on
end
semilogx(FreqFitabs_Hz, 20*log10(absRfit),'b')
hold off
grid on