%=============== settings.m ======================================
clear

addpath   ZZtoolbox/
addpath   ZZtoolbox/00pierrick/
% addpath ../00gabrielson/
%========= source of data
data_source = 'testbed_archive';
user        = 'charbit';
password    = 'sqlmomo';
channel     = '(''BDF'',''BDF'',''LWS'',''LWD'',''LKO'')';

yearstart   =  '2015';
monthstart  =  '06';
HMSstart    = '00:00:10';
yearend     =  '2015';
monthend    =  '06';
HMSend      = '23:50:10';

for ihc=1:1
    stations    = sprintf(' (''I26H%i'',''I26C%i'') ',ihc,ihc);
    for daystart_num    =  6 %1:2:25
        if daystart_num<10
            daystart    = ['0' num2str(daystart_num)];
            if daystart_num==9
                dayend  = '10';
            else
                dayend  = ['0' num2str(daystart_num+1)];
            end
        else
            daystart    = num2str(daystart_num);
            dayend      = num2str(daystart_num+1);
        end
        %=== extract data from the database
        h_starttime   = sprintf('%s/%s/%s %s',yearstart,monthstart,daystart, HMSstart);
        h_endtime     = sprintf('%s/%s/%s %s',yearend,monthend,dayend, HMSend);
        [~,starttime] = unix(['h2e ',h_starttime,' ofmt="%#"']);
        [~,endtime]   = unix(['h2e ',h_endtime,' ofmt="%#"']);
        starttime     = str2double(starttime);
        endtime       = str2double(endtime);
        wlength       = endtime-starttime;
        %========== temp fil in the current directory
        temporary_gparse_dir = 'ZZtoolbox/00pierrick/tempfiles/';
        %========== matlab format files are saved into the following directory
        savedirnamefull  ='../AAdataI26/';
        % %========= for analyzing data =================
        % MSCthreshold     = 0.97;
        % filtercharactfilename = 'filtercharacteristics1';
        pipeline
    end
end

%=================================================================

