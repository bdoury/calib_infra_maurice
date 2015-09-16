%=============== settings.m ======================================
addpath   ZZtoolbox/
% addpath ../5performCoherenceeffect/
% addpath ../00gabrielson/
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

%=== extract data from the database 
h_starttime   = sprintf('%s/%s/%s %s',yearstart,monthstart,daystart, HMSstart);
h_endtime     = sprintf('%s/%s/%s %s',yearend,monthend,dayend, HMSend);
[~,starttime] = unix(['h2e ',h_starttime,' ofmt="%#"']);
[~,endtime]   = unix(['h2e ',h_endtime,' ofmt="%#"']);
starttime     = str2double(starttime);
endtime       = str2double(endtime);
wlength       = endtime-starttime;

%========== temp fil in the current directory
temporary_gparse = 'tempfiles/gparse_temp.par';
filewfdisc       = 'tempfiles/gparse.wfdisc';
%========== matlab format files are saved into the following directory
savedirnamefull  ='/Users/maurice/etudes/ctbto/allJOBs2015/DATA_IS/I26/';
%========= for analyzing data =================
MSCthreshold     = 0.97;
filtercharactfilename = 'filtercharacteristics1';
%=================================================================

