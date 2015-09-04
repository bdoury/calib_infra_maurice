%============= pipeline.m =========================================
clear
%=== read settings 
settings
%=== extract data from the database 
h_starttime   = sprintf('%s/%s/%s %s',yearstart,monthstart,daystart, HMSstart);
h_endtime     = sprintf('%s/%s/%s %s',yearend,monthend,dayend, HMSend);
[~,starttime] = unix(['h2e ',h_starttime,' ofmt="%#"']);
[~,endtime]   = unix(['h2e ',h_endtime,' ofmt="%#"']);
starttime     = str2double(starttime);
endtime       = str2double(endtime);
wlength       = endtime-starttime;
fid           = fopen(temporary_gparse,'w');
%==== Write the query
fprintf(fid, 'open open data_source=%s user=%s password=%s\n', data_source, user, password);
fprintf(fid, '%s\n',['query wfdisc select * from sel3.wfdisc where sta in ', ...
    stations, 'and chan in ', channel,' and time between ',num2str(starttime),...
    ' and ',num2str(starttime+wlength),' order by sta,chan,time']);
fprintf(fid, '%s\n','read waveforms');
fprintf(fid, '%s\n','write waveforms');
fclose(fid);
%==== 
disp('***************************************************************')
disp('****************** query to data base *************************')
%==== Gparse run
unix('setenv ORACLE_HOME /cots/oracle/oracle-10.2;');
unix('setenv D_LIBRARY_PATH $ORACLE_HOME/lib:$ORACLE_HOME/lib32;');
unix('/ctbto/ims/sm/local/linux/Geotool++/2.3.10/bin/gparse<gparse_temp.par;');

%====
disp('***************************************************************')
disp('***************** convert to Matlab format ********************')
%====================== Convert to Matlab
filenamesavemat = convertCSStomatlab(filewfdisc,savedirnamefull);
%====
commandload     = sprintf('load %s',filenamesavemat);
eval(commandload)
%====
disp('***************************************************************')
disp('************************* analysis ****************************')
[SUTs, filteredsignals, allfrqsFFT_Hz, alltimes_sec, filterbank] = ...
    analyzeSUT(records,filtercharactfilename, MSCthreshold);
P = length(filterbank);
%==== to test display a figure
displayafigure
%====================== end ============================

