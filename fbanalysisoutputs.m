% Outputs:
%     SUTs: structures P x 1
%         xx.estimRsup: Rsup ratio structure
%         xx.estimRinf: Rinf ratio structure
%         xx.allMSCs: allMSCs structure
%         xx.Nsupthreshold: count of the number of values over
%                   the threshold
%         xx.Nsupthresholdintheband: counts of the number of values
%                   over the threshold in the filter bandwidth
%         xx.SCP = all spectral components
%         xx.indexinsidefreqband = P x 1, indices of the
%                    frequency bounds of each filter
%                    in the xx.frqsFFT_Hz
%         xx.theomodstdforRsup: theoretical STD of the module estimate
%         xx.stdPhase_rad: theoretical STD of the phase estimate in radian
%
%         alltimes_sec: cell  Px 1, each cell consists of
%             yy.FFT: time list in second of the DFTs
%             yy.SD: time list in second of the SCPs
%             yy.signals: time list in second of the signals
%
%         frqsFFT_Hz: cell P x 1, each cell consists of
%                   frequency list in Hz of the DFTs