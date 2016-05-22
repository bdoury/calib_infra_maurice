%====================== savesignals.m ====================================
% save data of a day pair into directorysave2daysignals
%==================== Warning ===============
% we have observed huge outliers in the following files:
% ihc==1, date==2015/10/07 ,from sample index 2.4e6
% ihc==1, date==2015/10/09 ,from sample index 2.4e6
% ihc==2, date==2015/08/07 ,from sample index 2.4e6
% ihc==2, date==2015/10/05 ,from sample index 2.4e6
% ihc==5, date==2015/10/05
% ihc==6, date==2015/10/07
% ihc==8, date==2015/10/07 , from sample index 2.5e6

% This program does not remove them for examination purposes
%
%============================================
setimesC_ihc                 = zeros(1,2);
setimesH_ihc                 = zeros(1,2);
problemHC                    = zeros(1,2);

idSc = 1;
idSh = 1;
cpC  = 1;
cpH  = 1;
idWS = 1;
idWD = 1;
idT  = 1;
signals     = zeros(34560000,2);
windSpeed   = zeros(34560000,1);
temperature = zeros(34560000,1);
windDir     = zeros(34560000,1);
Lrecords    = length(records);

for ir = 1:Lrecords
    switch records{ir}.channel
        case 'BDF'
            Fs_Hz = 20;%records{ir}.Fs_Hz;
            switch records{ir}.station(4)
                case 'C'
                    if cpC==1
                        cpC=cpC+1;
                        setimesC_ihc(1) = records{ir}.stime;
                    end
                    auxC = records{ir}.etime;
                    LLC = length(records{ir}.data);
                    signalsC = [records{ir}.data];
                    signals(idSc:idSc+LLC-1,2)=signalsC;
                    idSc = idSc+LLC;
                case 'H'
                    if cpH==1
                        cpH=cpH+1;
                        setimesH_ihc(1) = records{ir}.stime;
                    end
                    auxH = records{ir}.etime;
                    LLH = length(records{ir}.data);
                    signalsH = [records{ir}.data];
                    signals(idSh:idSh+LLH-1,1)=signalsH;
                    idSh = idSh+LLH;
            end
        case 'LKO'
            LLT = length(records{ir}.data);
            temperature(idT:idT+LLT-1) = [records{ir}.data];
            idT = idT + LLT;
        case 'LWS'
            LLWS = length(records{ir}.data);
            windSpeed(idWS:idWS+LLWS-1) = [records{ir}.data];
            idWS = idWS + LLWS;
        case 'LWD'
            FsWind_Hz = records{ir}.Fs_Hz;
            LLWD = length(records{ir}.data);
            Fs_wind_Hz = records{ir}.Fs_Hz;
            windDir(idWD:idWD+LLWD-1) = [records{ir}.data];
            idWD = idWD + LLWD;
    end
end
setimesC_ihc(2) = auxC;
setimesH_ihc(2) = auxH;

if not(idSc==idSh)
    fprintf('problem on %i\nLH = %i and LC = %i\n',...
        ihc,idSh,idSc)
    problemHC = [idSh,idSc];
end
idScMin = min([idSc, idSh]);

signals     = signals(1:idScMin-1,:);
windSpeed   = windSpeed(1:idWS-1);
windDir     = windDir(1:idWD-1);
temperature = temperature(1:idT-1);

Ts_sec = 1/Fs_Hz;
signals_centered = signals-ones(size(signals,1),1)*mean(signals);

commandsave = sprintf...
    ('save %ss%i/s%iyear%smonth%sday%s signals_centered Fs_Hz setimesC_ihc setimesH_ihc problemHC', ...
    directorysave2daysignals,ihc,ihc,yearstart,monthstart,daystart);
eval(commandsave)

if ihc==1
    commandsaveWIND = sprintf...
        ('save %ss%i/s%iWINDyear%smonth%sday%s windSpeed windDir temperature FsWind_Hz', ...
        directorysave2daysignals,ihc,ihc,yearstart,monthstart,daystart);
    eval(commandsaveWIND)
end



