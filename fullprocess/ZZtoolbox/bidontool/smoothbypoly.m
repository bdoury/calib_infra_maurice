function [HOut,rangefOut] = ...
    smoothbypoly(HIn,rangefIn,freqCut_Hz,LP)
%====================================================================
% The fit curve consists of 2 parts as polynoms.
% Inputs:
%       HIn : ordinates of the function
%       rangefIn : abscissas
%       freqCut_Hz : frequency where we change the polynom
% Outputs:
%       HOut : ordinates of the function
%       rangefOut : abscissas
%====================================================================

indexnan = isnan(HIn);
HIn = HIn(not(indexnan));
rangefOut = rangefIn(not(indexnan));
LrangefIn=length(rangefOut);
indexCut = find(rangefOut<freqCut_Hz,1,'last');
indexrangefreq1 = 1:indexCut;
indexrangefreq2 = indexCut:LrangefIn;
rangefreq1 = (rangefOut(indexrangefreq1))';
rangefreq2 = (rangefOut(indexrangefreq2))';
HH1 = exp(log(rangefreq1) * (0:LP(1)-1));
HH2 = exp(log(rangefreq2) * (0:LP(2)-1));
bet1 = HH1\HIn(indexrangefreq1);
bet2 = HH2\HIn(indexrangefreq2);
predHmedian1 = HH1*bet1;
predHmedian2 = HH2*bet2;

predHmedian1L = zeros(LrangefIn,1);
predHmedian2L = zeros(LrangefIn,1);
predHmedian1L(indexrangefreq1(1:indexCut-6))=predHmedian1(1:indexCut-6);
predHmedian2L(indexrangefreq2)=predHmedian2;
HOut = (predHmedian1L+predHmedian2L);
HOut(indexCut-6:indexCut) = ...
    ((6:-1:0)' .* predHmedian1(indexCut-6:indexCut)+(0:6)'*predHmedian2(1))/6;

