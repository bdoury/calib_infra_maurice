Pfilter=1;
filtercharact(Pfilter).designname     = 'butter';
filtercharact(Pfilter).Norder         = 4;
filtercharact(Pfilter).Wlow_Hz        = 0.8333;
filtercharact(Pfilter).Whigh_Hz       = 8;
filtercharact(Pfilter).SCPperiod_sec  = 30;
filtercharact(Pfilter).windowshape    = 'hann';
filtercharact(Pfilter).overlapDFT     = 0.5;
filtercharact(Pfilter).overlapSCP     = 0;
filtercharact(Pfilter).ratioDFT2SCP   = 5;
%------  
