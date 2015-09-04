%===========================================================
%         xx.Norder: order of the filter
%         xx.Wlow_Hz:
%         xx.Whigh_Hz:
%         xx.SCPperiod_sec: duration in second
%                 over which SPC is performed,
%                 expected as stationarity time
%         xx.windowshape: window shape for spectral
%                 analysis (ex. 'hann' or 'hamming', ...)
%         xx.overlapDFT: overlap rate for successive
%                 DFT (typically 0.5)
%         xx.overlapSCP: overlap rate for successive
%                 spectral components (typically 0)
%         xx.ratioDFT2SCP: ratio between period_sec and
%                 DFT duration (typical integer value is 5).
%===========================================================
P=1;
% filtercharact(P).designname = 'butter';
% filtercharact(P).Norder      = 2;
% filtercharact(P).Wlow_Hz       = 0.001;
% filtercharact(P).Whigh_Hz       = 0.2;
% filtercharact(P).SCPperiod_sec  = 2000;
% filtercharact(P).windowshape    = 'hann';
% filtercharact(P).overlapDFT     = 0.5;
% filtercharact(P).overlapSCP     = 0;
% filtercharact(P).ratioDFT2SCP  = 10;
% %------
% P=P+1;
% filtercharact(P).designname = 'butter';
% filtercharact(P).Norder      = 2;
% filtercharact(P).Wlow_Hz      = 0.2;
% filtercharact(P).Whigh_Hz      = 1;
% filtercharact(P).SCPperiod_sec      = 500;
% filtercharact(P).windowshape      = 'hann';
% filtercharact(P).overlapDFT      = 0.5;
% filtercharact(P).overlapSCP      = 0;
% filtercharact(P).ratioDFT2SCP     = 5;
% %------
% P=P+1;
filtercharact(P).designname = 'butter';
filtercharact(P).Norder      = 0;
filtercharact(P).Wlow_Hz      = 1;
filtercharact(P).Whigh_Hz      = 2;
filtercharact(P).SCPperiod_sec      = 1000;
filtercharact(P).windowshape      = 'hann';
filtercharact(P).overlapDFT      = 0.5;
filtercharact(P).overlapSCP      = 0;
filtercharact(P).ratioDFT2SCP     = 5;
% %------
% P=P+1;
% filtercharact(P).designname = 'butter';
% filtercharact(P).Norder      = 4;
% filtercharact(P).Wlow_Hz      = 2;
% filtercharact(P).Whigh_Hz      = 3;
% filtercharact(P).SCPperiod_sec      = 50;
% filtercharact(P).windowshape      = 'hann';
% filtercharact(P).overlapDFT      = 0.5;
% filtercharact(P).overlapSCP      = 0;
% filtercharact(P).ratioDFT2SCP     = 5;
% %------
% P=P+1;
% filtercharact(P).designname = 'butter';
% filtercharact(P).Norder      = 4;
% filtercharact(P).Wlow_Hz      = 3;
% filtercharact(P).Whigh_Hz      = 4;
% filtercharact(P).SCPperiod_sec       = 25;
% filtercharact(P).windowshape      = 'hann';
% filtercharact(P).overlapDFT      = 0.5;
% filtercharact(P).overlapSCP      = 0;
% filtercharact(P).ratioDFT2SCP     = 5;
% %------
% P=P+1;
% filtercharact(P).designname = 'butter';
% filtercharact(P).Norder      = 4;
% filtercharact(P).Wlow_Hz      = 4;
% filtercharact(P).Whigh_Hz      = 6;
% filtercharact(P).SCPperiod_sec     = 25;
% filtercharact(P).windowshape      = 'hann';
% filtercharact(P).overlapDFT      = 0.5;
% filtercharact(P).overlapSCP      = 0;
% filtercharact(P).ratioDFT2SCP     = 5;
% %------
% P=P+1;
% filtercharact(P).designname = 'butter';
% filtercharact(P).Norder      = 2;
% filtercharact(P).Wlow_Hz      = 3.5;
% filtercharact(P).Whigh_Hz      = 6;
% filtercharact(P).SCPperiod_sec     = 25;
% filtercharact(P).windowshape      = 'hann';
% filtercharact(P).overlapDFT      = 0.5;
% filtercharact(P).overlapSCP      = 0;
% filtercharact(P).ratioDFT2SCP      = 5;

% %------
% P=P+1;
% filtercharact(P).designname = 'butter';
% filtercharact(P).Norder      = 2;
% filtercharact(P).Wlow_Hz      = 5;
% filtercharact(P).Whigh_Hz      = 16;
% filtercharact(P).SCPperiod_sec     = 10;
% filtercharact(P).windowshape      = 'hann';
% filtercharact(P).overlapDFT      = 0.5;
% filtercharact(P).overlapSCP      = 0;
% filtercharact(P).ratioDFT2SCP      = 5;


% %------
% P=P+1;
% filtercharact(P).designname = 'butter';
% filtercharact(P).Norder      = 8;
% filtercharact(P).Wlow_Hz      = 4;
% filtercharact(P).Whigh_Hz      = 6;
% filtercharact(P).SCPperiod_sec     = 25;
% filtercharact(P).windowshape      = 'hann';
% filtercharact(P).overlapDFT      = 0.5;
% filtercharact(P).overlapSCP      = 0;
% filtercharact(P).ratioDFT2SCP      = 10;
