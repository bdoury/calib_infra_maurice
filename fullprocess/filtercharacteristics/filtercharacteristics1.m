Pfilter=1;
filtercharact(Pfilter).designname     = 'butter';
filtercharact(Pfilter).Norder         = 2;
filtercharact(Pfilter).Wlow_Hz        = 0.01;
filtercharact(Pfilter).Whigh_Hz       = 0.06;
filtercharact(Pfilter).SCPperiod_sec  = 1000;
filtercharact(Pfilter).windowshape    = 'hann';
filtercharact(Pfilter).overlapDFT     = 0.5;
filtercharact(Pfilter).overlapSCP     = 0;
filtercharact(Pfilter).ratioDFT2SCP   = 5;
%------  
Pfilter=Pfilter+1;
filtercharact(Pfilter).designname     = 'butter';
filtercharact(Pfilter).Norder         = 2;
filtercharact(Pfilter).Wlow_Hz        = 0.05;
filtercharact(Pfilter).Whigh_Hz       = 0.15;
filtercharact(Pfilter).SCPperiod_sec  = 700;
filtercharact(Pfilter).windowshape    = 'hann';
filtercharact(Pfilter).overlapDFT     = 0.5;
filtercharact(Pfilter).overlapSCP     = 0;
filtercharact(Pfilter).ratioDFT2SCP   = 5;
%------  
Pfilter=Pfilter+1;
filtercharact(Pfilter).designname     = 'butter';
filtercharact(Pfilter).Norder         = 2;
filtercharact(Pfilter).Wlow_Hz        = 0.14;
filtercharact(Pfilter).Whigh_Hz       = 0.35;
filtercharact(Pfilter).SCPperiod_sec  = 250;
filtercharact(Pfilter).windowshape    = 'hann';
filtercharact(Pfilter).overlapDFT     = 0.5;
filtercharact(Pfilter).overlapSCP     = 0;
filtercharact(Pfilter).ratioDFT2SCP   = 5;
%------  
Pfilter=Pfilter+1;
filtercharact(Pfilter).designname     = 'butter';
filtercharact(Pfilter).Norder         = 2;
filtercharact(Pfilter).Wlow_Hz        = 0.32;
filtercharact(Pfilter).Whigh_Hz       = 0.7;
filtercharact(Pfilter).SCPperiod_sec  = 150;
filtercharact(Pfilter).windowshape    = 'hann';
filtercharact(Pfilter).overlapDFT     = 0.5;
filtercharact(Pfilter).overlapSCP     = 0;
filtercharact(Pfilter).ratioDFT2SCP   = 5;
%------  
Pfilter=Pfilter+1;
filtercharact(Pfilter).designname     = 'butter';
filtercharact(Pfilter).Norder         = 3;
filtercharact(Pfilter).Wlow_Hz        = 0.6;
filtercharact(Pfilter).Whigh_Hz       = 1.5;
filtercharact(Pfilter).SCPperiod_sec  = 30;
filtercharact(Pfilter).windowshape    = 'hann';
filtercharact(Pfilter).overlapDFT     = 0.5;
filtercharact(Pfilter).overlapSCP     = 0;
filtercharact(Pfilter).ratioDFT2SCP   = 5;
%------  
Pfilter=Pfilter+1;
filtercharact(Pfilter).designname     = 'butter';
filtercharact(Pfilter).Norder         = 4;
filtercharact(Pfilter).Wlow_Hz        = 1.5;
filtercharact(Pfilter).Whigh_Hz       = 6;
filtercharact(Pfilter).SCPperiod_sec  = 20;
filtercharact(Pfilter).windowshape    = 'hann';
filtercharact(Pfilter).overlapDFT     = 0.5;
filtercharact(Pfilter).overlapSCP     = 0;
filtercharact(Pfilter).ratioDFT2SCP   = 5;
%------  
