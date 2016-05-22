function [BUT,AUT] = extractHUTols(SRR,SUU,SUR,Q,P, BREF, AREF, flagbadUU)
%================================================
% 1 = REF and UT = 2
Lfft = size(SRR,1);
if flagbadUU
    SB = AREF .* SRR;
    SA = BREF .* SUR;
    FB = fft(eye(Lfft,Q),Lfft);
    FA = fft(eye(Lfft,P),Lfft);
    TB = FB .* (SB*ones(1,Q));
    TA = FA(:,2:P) .* (SA*ones(1,P-1));
    T = [TB -TA];
    beat = T\SA;
    BUT = real(beat(1:Q));
    AUT = real([1 ; beat(Q+1:P+Q-1)]);
else
    SA = SUU .* BREF;
    SB = SUR .* AREF;
    FB = fft(eye(Lfft,Q),Lfft);
    FA = fft(eye(Lfft,P),Lfft);
    TB = FB .* (SB*ones(1,Q));
    TA = FA(:,2:P) .* (SA*ones(1,P-1));
    T = [TB -TA];
    beat = (T)\SA;
    BUT = real(beat(1:Q));
    AUT = real([1 ; beat(Q+1:P+Q-1)]);
end
%=========



