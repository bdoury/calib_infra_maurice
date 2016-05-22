clear
%==========================================================
% this programs calls
%   . statsonH11H12.m
%   . allprintstatsonHest
% Perform the RMSE as MSC based on the spectral matrix
% Theoretical distributions are used.
% No signal is generated. 
% For comparison the progs estimHanalysis.m generated 
% signals to do the same calculation.
%==========================================================
absHUonHR              = 1;
statsonH11H12
cdesave = sprintf('save statsonHest%3.2f.mat',absHUonHR);
eval(cdesave)
% allprintstatsonHest
%==========================================================