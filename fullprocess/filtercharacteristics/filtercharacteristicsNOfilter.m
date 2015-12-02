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
Pfilter=1;
filtercharact(Pfilter).designname     = 'butter';
filtercharact(Pfilter).Norder         = 0;
filtercharact(Pfilter).Wlow_Hz        = 0.001;
filtercharact(Pfilter).Whigh_Hz       = 10;
filtercharact(Pfilter).SCPperiod_sec  = 500;
filtercharact(Pfilter).windowshape    = 'hann';
filtercharact(Pfilter).overlapDFT     = 0.5;
filtercharact(Pfilter).overlapSCP     = 0;
filtercharact(Pfilter).ratioDFT2SCP   = 5;
%------