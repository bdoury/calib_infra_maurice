function [swelch,s] = mypwelch(x,y,Lfft,overlap)


N = length(x);
tra = fix(Lfft*(1-overlap));
nblocks = fix(N/tra);
x=[x;zeros(Lfft,1)];
y=[y;zeros(Lfft,1)];
s=zeros(Lfft,nblocks);
windowHam = hamming(Lfft);
for in = 0:nblocks-1
    it1=in*tra+1;it2=it1+Lfft-1;
    x1=x(it1:it2);
    y1=y(it1:it2);
    x1=x1 .* windowHam;
    y1=y1 .* windowHam;
    sx1=fft(x1,Lfft);
    sy1=fft(y1,Lfft);
    s(:,in+1)=sx1 .* conj(sy1);
end
swelch = mean(s,2)/Lfft;
%======================

    