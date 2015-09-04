function [Xp, times_blocks] = ...
    computeXfMat(xREF,xUT,Lfft, overrate, Fs)

x=[xREF xUT];
[lx,nbsignals] = size(x);
trans = fix(Lfft*(1-overrate));
nbblocks = fix(lx/trans);
times_blocks.pts = ((0:nbblocks-1)+1/2);
times_blocks.sec = times_blocks.pts*trans/Fs;

Xp = zeros(Lfft,nbsignals,nbblocks);
for ib = 0:nbblocks-2
    it1  = ib*trans+1;
    it2  = min([it1+Lfft-1, lx]);
    UU   = ones(it2-it1+1,1);
    xaux = x(it1:it2,:);
    xaux = xaux-UU*mean(xaux,1);
    Wwind = hamming(it2-it1+1,'periodic');    
    for i1=1:nbsignals
        xwinaux = xaux(:,i1) .* Wwind;
        Xp(:,i1,ib+1) = fft(xwinaux,Lfft)/sqrt(Lfft);
    end
end