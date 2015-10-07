clear

id6=SUTs(6).allMSCs.indexcst(idipinf(6):idipsup(6),:);
 
  ind6=find(sum(id6,1)==4)
 selind = 972;
 T6 = alltimes_sec{6}.SD(selind);
n6=T6*20;
xUT6=filteredsignals(n6-1000:n6+1000,2,6);
xRF6=filteredsignals(n6-1000:n6+1000,1,6);
plot(alltimes_sec{6}.signals(n6-1000:n6+1000),[xUT6 xRF6])
20*log10(std(xUT6) / std(xRF6))

[PUT,W] = pwelch(xUT6,[],[],[],Fs_Hz);
[PRF,W] = pwelch(xRF6,[],[],[],Fs_Hz);

%  SUTs(6).allMSCs.tab(idipinf(6):idipsup(6),selind)

plot(W,10*log10(PUT))
hold on
plot(W,10*log10(PRF),'r')
hold off