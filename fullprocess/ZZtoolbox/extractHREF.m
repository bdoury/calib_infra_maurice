function Rref2ut = extractHREF(xUT, xREF, Lfft, g, overlapFFT)


T = length(xREF);
shiftFFT = (1-overlapFFT)*Lfft;
nbFrames = fix((T-Lfft)/shiftFFT);
allFFTsUT = zeros(Lfft,nbFrames);
allFFTsREF = zeros(Lfft,nbFrames);
wh = hann(Lfft,'periodic');
for ib=0:nbFrames-1
    id1=ib*shiftFFT+1;
    id2=id1+Lfft-1;
    id2 = min([id2, T]);
    allFFTsREF(:,ib+1)=fft(xREF(id1:id2) .* wh)/sqrt(Lfft);
    allFFTsUT(:,ib+1)=fft(xUT(id1:id2) .* wh)/sqrt(Lfft);
end
meanFFTUU = mean(abs(allFFTsUT).^2,2);
meanFFTRR = mean(abs(allFFTsREF).^2,2);
meanFFTUR = mean((allFFTsUT .* conj(allFFTsREF)),2);
for ib=1:Lfft-1
    plot(allFFTsREF(1+ib,:),'.')
end

gm1 = g-1;
Rref2ut = (-gm1* (abs(meanFFTUR)) + ...
    sqrt(gm1*gm1*(abs(meanFFTUR)) .^2 ...
    +4*g*meanFFTUU.*meanFFTRR)) ./ (2*meanFFTRR);
    