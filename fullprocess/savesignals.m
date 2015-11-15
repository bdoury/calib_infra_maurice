%====================== savesignals.m ====================================
% the data are read form a file with the name built as:
%                  sta1_Y2015_D239.mat
% meaning station number 1, year 2015, days 239 and 240
% but the data are in IDC format as
% They  have been extracted from IDC testbed for 2 days.
% The signals are saved into the directory directorysave2daysignals
% for wind in station data are saved.
%=========================================================================
clear
allcolors = ['b.';'r.';'m.';'c.';'g.';'k.';'rx';'yx';'mx';'rx';'kx';...
    'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];

%=============== Warning ===============
% we have observed huge outliers in the following files:
% ihc==2, ifile== 33, i.e. sta2_Y2015_D219.mat from sample index 2.4e6
% ihc==8, ifile==63, i.e. sta8_Y2015_D280.mat from sample index 2.5e6
% ihc==5,ifile==62
% ihc==6,ifile==63
% This program does not remove them for examination purposes 
%
%============================================
%============================================
%============================================
%=====================
MSCthreshold = 0.98;
%=====================
ihc   = 3;
directorydatafromIDC  = '../../../AAdataI26formatIDC/';
% directory to save data
directorysave2daysignals  = '../../../AAdataI26calib/';
%===================== read data =========================
fileswithdotmat              = dir(sprintf('%ss%i/sta%i*.mat', ...
    directorydatafromIDC,ihc,ihc));
nbmats                       = length(fileswithdotmat);

setimesC_ihc                 = zeros(1,2);
setimesH_ihc                 = zeros(1,2);
problemHC                    = zeros(1,2);
%==================================================
for ifile = 1:nbmats,ifile
    fullfilename_i      = fileswithdotmat(ifile).name;
    dotlocation         = strfind(fullfilename_i,'.');
    underscorelocation  = strfind(fullfilename_i,'_');
    filenameonly        = fullfilename_i(setdiff(1:dotlocation-1,...
        underscorelocation));
    commandload         = sprintf('load %ss%i/%s',...
        directorydatafromIDC,ihc,fullfilename_i);
    eval(commandload)
    
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
                            setimesC_ihc(ifile,1) = records{ir}.stime;
                        end
                        auxC = records{ir}.etime;
                        LLC = length(records{ir}.data);
                        signalsC = [records{ir}.data];
                        signals(idSc:idSc+LLC-1,2)=signalsC;
                        idSc = idSc+LLC;
                    case 'H'
                        if cpH==1
                            cpH=cpH+1;
                            setimesH_ihc(ifile,1) = records{ir}.stime;
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
    setimesC_ihc(ifile,2) = auxC;
    setimesH_ihc(ifile,2) = auxH;
    
    if not(idSc==idSh)
        fprintf('problem with file # %i on %i\nLH = %i and LC = %i\n',...
            ifile,ihc,idSh,idSc)
        problemHC(ifile,:) = [idSh,idSc];
    end
    idScMin = min([idSc, idSh]);
    
    signals     = signals(1:idScMin-1,:);
    
    
    windSpeed   = windSpeed(1:idWS-1);
    windDir     = windDir(1:idWD-1);
    temperature = temperature(1:idT-1);
    
    Ts_sec = 1/Fs_Hz;
    signals_centered = signals-ones(size(signals,1),1)*mean(signals);
    
    
    commandsave = sprintf...
  ('save %ss%i/year%sday%s signals_centered Fs_Hz setimesC_ihc problemHC', ...
        directorysave2daysignals,ihc,...
        filenameonly(6:9),filenameonly(11:13));
    eval(commandsave)
    
    if ihc==1
        commandsaveWIND = sprintf...
 ('save %ss%i/WINDyear%sday%s windSpeed windDir temperature FsWind_Hz', ...
        directorysave2daysignals,ihc,...
        filenameonly(6:9),filenameonly(11:13));
    eval(commandsaveWIND)
    end  
end



