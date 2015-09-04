function Sxyf = crosssp(x1, x2, Lfft, overlap)
%===================================================
%
%
%
%
%===================================================
x1=x1(:); x2=x2(:);
N = length(x1);
shiftFFT = fix(Lfft*(1-overlap))+1;
nbblocksFFT = fix(N/shiftFFT);
Sxyf = zeros(Lfft,1);
hanwin = hann(Lfft,'periodic');
for ib=0:nbblocksFFT-1,
    it1=ib*shiftFFT+1;
    it2=it1+Lfft-1;
    if it2>N
        it2 = N;
        x1e = [x1(it1:it2);zeros(Lfft-(it2-it1)-1,1)];
        x2e = [x2(it1:it2);zeros(Lfft-(it2-it1)-1,1)];
    else
        x1e = x1(it1:it2);
        x2e = x2(it1:it2);
    end
    x1wf = fft(x1e .* hanwin,Lfft);
    x2wf = fft(x2e .* hanwin,Lfft);
    Sxyf = Sxyf+x1wf .* conj(x2wf);
end
Sxyf = Sxyf/nbblocksFFT/sqrt(Lfft);
%===================================================