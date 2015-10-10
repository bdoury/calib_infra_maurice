%====================== estimationwithFB.m =============================
clear
allcolors = ['b.';'r.';'m.';'c.';'g.';'k.';'rx';'yx';'mx';'rx';'kx';...
    'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];

addpath ZZtoolbox/

directoryresults      = 'AAresultswithFBXsepar';
directorydatafromIDC  = '../../../AAdataI26/';

%==========================================================================
% the data are in a file with the name built as:
%                  sta1_Y2015_D239.mat
% meaning station number 1, year 2015, days 239 and 240
% these data have been extracted from IDC testbed
%==========================================================================

%============== run the filter bank characteristics

%============== run the filter bank characteristics
% OLD execute
% filtercharactfilename = 'filtercharacteristics.m';
%======
% generate by geneFB.m, use LOAD
filtercharactfilename = 'filtercharacteristics2.mat';
cmdloadcharact = sprintf('load(''%s'')',filtercharactfilename);

%======
% to select only 1Hz
% filtercharactfilename = 'filtercharacteristics0812Hz.m';
% cmdloadcharact = sprintf('run(''%s'')',filtercharactfilename);

%======
eval(cmdloadcharact);
Pfilter = length(filtercharact);
if and(Pfilter==1, filtercharact(Pfilter).Norder==0)
    filtercharact(Pfilter).Wlow_Hz  = 0.001;
    filtercharact(Pfilter).Whigh_Hz = 10;
end
%=====================
MSCthreshold   = 0.98;
%=====================
for indexofSTA = 4
    %=====================
    % under test = 1, reference = 2
    %===================== read data =========================
    fileswithdotmat           = dir(sprintf('%ss%i/sta%i*.mat',directorydatafromIDC,indexofSTA,indexofSTA));
    nbmats                    = length(fileswithdotmat);
    
    for ifile=1:nbmats, ifile, tic
        
        allRatioPfilters          = zeros(10000,1);
        allfrqsPfilters           = zeros(10000,1);
        allSTDmodRatioPfilters    = zeros(10000,1);
        allSTDphaseRatioPfilters  = zeros(10000,1);
        allmeanMSCcstPfilters     = zeros(10000,1);
        nbofvaluesoverthreshold   = zeros(10000,1);

        fullfilename_i      = fileswithdotmat(ifile).name;
        dotlocation         = strfind(fullfilename_i,'.');
        underscorelocation  = strfind(fullfilename_i,'_');
        filenameonly        = fullfilename_i(setdiff(1:dotlocation-1,underscorelocation));
        commandload         = sprintf('load %ss%i/%s',directorydatafromIDC,indexofSTA,fullfilename_i);
        eval(commandload)
        
        idSc = 1;
        idSh = 1;
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
        size(signals)
        windSpeed = windSpeed(1:idWS-1);
        windDir = windDir(1:idWD-1);
        temperature = temperature(1:idT-1);
        
        Ts_sec = 1/Fs_Hz;
        signals_centered=signals-ones(size(signals,1),1)*mean(signals);
        
        %============================================
        % notice that the SUTs is not saved, therefore we have only the
        % last associated the laxt index which is NBMATS
        [SUTs, filteredsignals, allfrqsFFT_Hz, alltimes_sec, filterbank] = ...
            fbankanalysis(signals_centered,filtercharact,Fs_Hz,MSCthreshold);
        
        %============================================
        P       = length(SUTs);
        idipinf = zeros(P,1);
        idipsup = zeros(P,1);
        id1     = 1;
        for ip=1:P
            idipinf(ip) = SUTs(ip).indexinsidefreqband(1);
            idipsup(ip) = SUTs(ip).indexinsidefreqband(2);
            id2         = id1+(idipsup(ip)-idipinf(ip));
            allRatioPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRsup.modcst(idipinf(ip):idipsup(ip)) .* ...
                exp(1i*SUTs(ip).estimRsup.phasecst(idipinf(ip):idipsup(ip)));
            allSTDmodRatioPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRsup.stdmodcst(idipinf(ip):idipsup(ip));
            allSTDphaseRatioPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRsup.phasecst(idipinf(ip):idipsup(ip));
            allmeanMSCcstPfilters(id1:id2,ifile) = ...
                nanmean(SUTs(ip).allMSCs.tabcst(idipinf(ip):idipsup(ip),:),2);
            allfrqsPfilters(id1:id2,ifile) = ...
                SUTs(ip).frqsFFT_Hz(idipinf(ip):idipsup(ip))';
            nbofvaluesoverthreshold(id1:id2,ifile) = ...
                sum(not(isnan(SUTs(ip).allMSCs.tabcst(idipinf(ip):idipsup(ip),:))),2);
            id1 = id2+1;
        end
        allRatioPfilters         = allRatioPfilters(1:id1-1,:);
        allfrqsPfilters          = allfrqsPfilters(1:id1-1,1);
        allSTDmodRatioPfilters   = allSTDmodRatioPfilters(1:id1-1,:);
        allSTDphaseRatioPfilters = allSTDphaseRatioPfilters(1:id1-1,:);
        allmeanMSCcstPfilters    = allmeanMSCcstPfilters(1:id1-1,:);
        nbofvaluesoverthreshold  = nbofvaluesoverthreshold(1:id1-1,:);
        
        
        Dstart = str2double(fileswithdotmat(ifile).name(13:15));

        comsave          = ...
            sprintf('save %s/s%i/resultssta26sensor%iday%i',directoryresults,indexofSTA,indexofSTA,Dstart);
        clear signals
        clear signalsC
        clear signalsH
        clear signals_centered
        clear filteredsignals
        clear windDir
        clear windSpeed
        clear records
        clear temperature
        clear alltimes_sec
        clear SUTs
        eval(comsave);
    end
end
%============================ END =========================================