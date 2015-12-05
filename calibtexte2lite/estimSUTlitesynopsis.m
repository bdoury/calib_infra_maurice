function [Rsup, freqslin, STDmoduleRlin, ...
    STDphaseRlin_rd, nboverTHlin] = ...
    estimSUTlite ...
    (signals, structfiltercharacteristics, frequencylist_Hz, ...
    Fs_Hz, MSCthreshold, trimpercent)
%========================================================================
% Synopsis:
% [Rsup, freqslin, STDmoduleRlin, ...
%     STDphaseRlin_rd, nboverTHlin] = ...
%     estimSUTlite ...
%     (signals, structfiltercharacteristics, frequencylist_Hz, ...
%     Fs_Hz, MSCthreshold, trimpercent)
%===============
% Inputs:
%     - signals : T x 2
%     - structfiltercharacteristics (FB structure)
%           see document
%     - frequencylist_Hz: array N x 1 of the selected frequencies
%       in Hz. N can take any value under Fs_Hz/2
%     - Fs_Hz: sampling frequency in Hz
%     - MSCthreshold:
%     - trimpercent: percent of values keptfor averaging
%===============
% Outputs:
%     - Rsup: array N x 1 of the estimated ratios
%     - freqslin: array N x 1 of the selected frequencies
%       in Hz. almost the same as frequencylist_Hz, except if some 
%       are outside of the FB bandwidths.
%     - STDmoduleR: array N x 1 of the STD on the module of Rsup
%     - STDphaseR_rd: array N x 1 of the STD on the phase of Rsup
%     - nboverTH: array N x 1 of the number of values over the threshold
%========================================================================