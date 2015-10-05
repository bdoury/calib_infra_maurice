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
filtercharact(Pfilter).Norder         = 2;
filtercharact(Pfilter).Wlow_Hz        = 0.005;
filtercharact(Pfilter).Whigh_Hz       = 0.0175;
filtercharact(Pfilter).SCPperiod_sec  = 1000;
filtercharact(Pfilter).windowshape    = 'hann';
filtercharact(Pfilter).overlapDFT     = 0.5;
filtercharact(Pfilter).overlapSCP     = 0;
filtercharact(Pfilter).ratioDFT2SCP   = 5;
%------
Pfilter=Pfilter+1;
filtercharact(Pfilter).designname       = 'butter';
filtercharact(Pfilter).Norder           = 3;
filtercharact(Pfilter).Wlow_Hz          = 0.0175;
filtercharact(Pfilter).Whigh_Hz         = 0.0613;
filtercharact(Pfilter).SCPperiod_sec    = 300;
filtercharact(Pfilter).windowshape      = 'hann';
filtercharact(Pfilter).overlapDFT       = 0.5;
filtercharact(Pfilter).overlapSCP       = 0;
filtercharact(Pfilter).ratioDFT2SCP     = 5;
%------
Pfilter=Pfilter+1;
filtercharact(Pfilter).designname       = 'butter';
filtercharact(Pfilter).Norder           = 3;
filtercharact(Pfilter).Wlow_Hz          = 0.0613;
filtercharact(Pfilter).Whigh_Hz         = 0.2144;
filtercharact(Pfilter).SCPperiod_sec    = 80;
filtercharact(Pfilter).windowshape      = 'hann';
filtercharact(Pfilter).overlapDFT       = 0.5;
filtercharact(Pfilter).overlapSCP       = 0;
filtercharact(Pfilter).ratioDFT2SCP     = 5;
%------
Pfilter=Pfilter+1;
filtercharact(Pfilter).designname       = 'butter';
filtercharact(Pfilter).Norder           = 3;
filtercharact(Pfilter).Wlow_Hz          = 0.2144;
filtercharact(Pfilter).Whigh_Hz         = 0.75;
filtercharact(Pfilter).SCPperiod_sec    = 25;
filtercharact(Pfilter).windowshape      = 'hann';
filtercharact(Pfilter).overlapDFT       = 0.5;
filtercharact(Pfilter).overlapSCP       = 0;
filtercharact(Pfilter).ratioDFT2SCP     = 5;
%------
Pfilter=Pfilter+1;
filtercharact(Pfilter).designname       = 'butter';
filtercharact(Pfilter).Norder           = 3;
filtercharact(Pfilter).Wlow_Hz          = 0.75;
filtercharact(Pfilter).Whigh_Hz         = 3;
filtercharact(Pfilter).SCPperiod_sec    = 10;
filtercharact(Pfilter).windowshape      = 'hann';
filtercharact(Pfilter).overlapDFT       = 0.5;
filtercharact(Pfilter).overlapSCP       = 0;
filtercharact(Pfilter).ratioDFT2SCP     = 5;
%------
Pfilter=Pfilter+1;
filtercharact(Pfilter).designname       = 'butter';
filtercharact(Pfilter).Norder           = 3;
filtercharact(Pfilter).Wlow_Hz          = 3;
filtercharact(Pfilter).Whigh_Hz         = 9;
filtercharact(Pfilter).SCPperiod_sec    = 5;
filtercharact(Pfilter).windowshape      = 'hann';
filtercharact(Pfilter).overlapDFT       = 0.5;
filtercharact(Pfilter).overlapSCP       = 0;
filtercharact(Pfilter).ratioDFT2SCP     = 5;
%------

