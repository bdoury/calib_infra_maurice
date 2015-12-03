clear 
overlapSD=0.5;

Fs_Hz = 20;
N = 24*3600*20;
TSCP_sec = 1000;
NaverageFFTs = 5;

Lfft = TSCP_sec*Fs_Hz/NaverageFFTs;

overlapFFT = 0.5; 


shiftSignal = fix((1-overlapFFT)*Lfft);
NblocksFFT  = fix((N-(Lfft-shiftSignal))/shiftSignal);

NaverageFFT = fix(NaverageFFTs/(1-overlapFFT))-1;
shiftFFTs = fix((1-overlapSD)*NaverageFFTs);
NshiftFFTs_with_overlap = max([shiftFFTs/(1-overlapFFT)-1,1]);
NSD       = fix(NblocksFFT/NshiftFFTs_with_overlap);

for ibB=1:NSD,
    indB1 = (ibB-1)*NshiftFFTs_with_overlap+1;
    indB2 = indB1+NaverageFFT-1;
    indB  = fix(indB1):fix(indB2);
    indBsave{ibB}  = indB(indB<= NblocksFFT);
end
