%=============== settings.m ======================================
addpath  ../../../ZZtoolbox/
addpath ../5performCoherenceeffect/
addpath ../00gabrielson/
%========= source of data
data_source = 'testbed_archive';
user        = 'charbit';
password    = 'sqlmomo';
channel     = '(''BDF'',''BDF'',''LWS'',''LWD'',''LKO'')';
stations    = '(''I26H6'',''I26C6'')';
yearstart   =  '2015';
monthstart  =  '06';
daystart    =  '06';
HMSstart    = '00:00:10';
yearend     =  '2015';
monthend    =  '06'; 
dayend      =  '08';
HMSend      = '23:50:10';

%========== temp fil in the current directory
temporary_gparse = 'gparse_temp.par';
filewfdisc       = 'gparse.wfdisc';
%========== matlab format files are saved into the following directory
savedirnamefull  ='/dvlscratch/SHI/users/charbit/ProjectIMS2015b/DATA_IS/I26/';
%========= for analyzing data =================
MSCthreshold     = 0.97;
filtercharactfilename = 'filtercharacteristics';
%=================================================================

