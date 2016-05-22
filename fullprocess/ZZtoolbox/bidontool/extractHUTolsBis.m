function [BUT,AUT] = extractHUTolsBis(Cww,Cwx,Q,P,W)
%================================================
% 1 = REF and UT = 2
if nargin == 4, W = 1;end
Lfft = size(Cww,1);
FB = fft(eye(Lfft,Q),Lfft);
FA = fft(eye(Lfft,P),Lfft);
TB = FB .* (Cww*ones(1,Q));
TA = FA .* (Cwx*ones(1,P));
T = [TB -TA(:,2:P)];
beat = (W*T)\(W*TA(:,1));
BUT = real(beat(1:Q));
AUT = real([1 ; beat(Q+1:P+Q-1)]);

%=========



