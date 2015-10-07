%====================== extractSensorRatio.m ================
% we only compute the response at one frequency
% determined by the filter, typically 1 Hz
%======================
clear
allcolors = ['b.';'r.';'m.';'c.';'g.';'k.';'rx';'yx';'mx';'rx';'kx';...
    'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];
%=====================
%


directoryresults = 'AAresults0812Hz';
%
%=========================================================================
addpath ZZtoolbox/00gabrielson/
addpath ZZtoolbox/00benoit/
%=======
%
directorydata = '../../../AAdataI26/';
%
%=========================================================================
addpath ZZtoolbox/

%======
% to select only 1Hz
filtercharactfilename = 'filtercharacteristics0812Hz.m';
cmdloadcharact = sprintf('run(''%s'')',filtercharactfilename);
eval(cmdloadcharact);

%=====================
MSCthreshold = 0.98;
%=====================
for indexofSTA = 1:8,indexofSTA
    %=====================
    % under test = 1, reference = 2
    %===================== read data =========================
    filesmat = dir(sprintf('%ss%i/sta%i*.mat',directorydata,indexofSTA,indexofSTA));
    nbmats   = length(filesmat);
    allRatioPfilters          = zeros(10000,nbmats);
    allfrqsPfilters           = zeros(10000,nbmats);
    allSTDmodRatioPfilters    = zeros(10000,nbmats);
    allSTDphaseRatioPfilters  = zeros(10000,nbmats);
    allmeanMSCcstPfilters     = zeros(10000,nbmats);
    nbofvaluesoverthreshold   = zeros(10000,nbmats);
    
    for ifile=1:nbmats,ifile
        file_i=filesmat(ifile).name;
        dotpos = strfind(file_i,'.');
        underscorepos = strfind(file_i,'_');
        fileonly = file_i(setdiff(1:dotpos-1,underscorepos));
        commandload    = sprintf('load %ss%i/%s',directorydata,indexofSTA,file_i);
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
                
        Lrecords =length(records);
        for ir =1:Lrecords
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
        windSpeed = windSpeed(1:idWS-1);
        windDir = windDir(1:idWD-1);
        temperature = temperature(1:idT-1);
        
        Ts_sec = 1/Fs_Hz;
        signals_centered=signals-ones(size(signals,1),1)*mean(signals);
        
        %============================================
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
            id2 = id1+(idipsup(ip)-idipinf(ip));
            allRatioPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRsup.modcst(idipinf(ip):idipsup(ip)) .* ...
                exp(1i*SUTs(ip).estimRsup.phasecst(idipinf(ip):idipsup(ip)));
            allfrqsPfilters(id1:id2,ifile) = ...
                SUTs(ip).frqsFFT_Hz(idipinf(ip):idipsup(ip))';
            id1 = id2+1;
        end
        
    end
    allRatioPfilters         = allRatioPfilters(1:id1-1,:);
    allfrqsPfilters          = allfrqsPfilters(1:id1-1,1);
    allSTDmodRatioPfilters   = allSTDmodRatioPfilters(1:id1-1,:);
    allSTDphaseRatioPfilters = allSTDphaseRatioPfilters(1:id1-1,:);
    allmeanMSCcstPfilters    = allmeanMSCcstPfilters(1:id1-1,:);
    comsave          = ...
        sprintf('save %s/resultssta26sensor%i%s',directoryresults,indexofSTA);
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
    allRatioPfilters = allRatioPfilters(1:id1-1,:);
    allfrqsPfilters  = allfrqsPfilters(1:id1-1,1);
    comsave          = ...
        sprintf('save %s/resultssta26sensor%i%s',directoryresults,indexofSTA);
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
%============================ END =========================================