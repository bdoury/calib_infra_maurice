function [SUTs, filteredsignals, allfrqsFFT_Hz, ...
    alltimes_sec, filterbank] = ...
    analyzeSUT(records,filtercharactfilename, MSCthreshold)
%=================================================================
% This function uses the "records" extracted from the database
% by the program "convertCSStomatlab.m".
<<<<<<< HEAD:fullprocess/ZZtoolbox/analyzeSUT.m
% This function then calls the function "fbankanalysis.m" 
=======
% This function thencalls the function "fbankanalysis.m" 
>>>>>>> e518a59897f74aa44d579d6b4e37dafce4bba1aa:fullprocess/analyzeSUT.m
% of the ZZtoolbox which returns the structure "SUTs".
% This structure is of general purposes. The program
% "displayafigure.m" is an example.
%=================================================================

idSc = 1;
idSh = 1;
idWS = 1;
idWD = 1;
idT  = 1;
signals     = zeros(34560000,2);
windSpeed   = zeros(34560000,1);
temperature = zeros(34560000,1);
windDir     = zeros(34560000,1);

Lrecords =length(records);
for ir =1:Lrecords
    switch records{ir}.channel
        case 'BDF'
            Fs_Hz = records{ir}.Fs_Hz;
            switch records{ir}.station(4)
                case 'C'
                    LLC = length(records{ir}.data);
                    signalsC = [records{ir}.data];
                    signals(idSc:idSc+LLC-1,2)=signalsC;
                    idSc = idSc+LLC;
                case 'H'
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
signals = signals(1:idSc-1,:);
windSpeed = windSpeed(1:idWS-1);
windDir = windDir(1:idWD-1);
temperature = temperature(1:idT-1);
Ts_sec = 1/Fs_Hz;
signals_centered=signals-ones(size(signals,1),1)*mean(signals);
%============================================
cmdloadcharact = sprintf('run(''%s'')',filtercharactfilename);
eval(cmdloadcharact);
%============================================
[SUTs, filteredsignals, allfrqsFFT_Hz, alltimes_sec, filterbank] = ...
    fbankanalysis(signals_centered,filtercharact,Fs_Hz,MSCthreshold);
%============================================