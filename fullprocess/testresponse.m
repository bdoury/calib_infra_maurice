clear all
close all

addpath ZZtoolbox/00benoit/

freq_vector = logspace(log10(0.001),log10(30),300) .';


TFMB2005     = idc2fap_is(freq_vector,'I26DE_BDF_RSP_2015134_MB2005');
TFMB3        = idc2fap_is(freq_vector,'I26DE_BDF_RSP_2015134_MB3');

plot(freq_vector,(abs([TFMB3 TFMB2005])));
set(gca, 'ylim',[0.9 1.2])
set(gca, 'xscale','log')
grid on
