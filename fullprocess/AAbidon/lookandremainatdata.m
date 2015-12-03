% function remainindex = lookandremainatdata...
%     (directorydatafromIDC,ihc, flaglook)
%===============================================================

% directorydatafromIDC: directory where are saved the
% signals in .mat. This directiory contains 8 directories
% from s1 to s8.
% ihc: sensor number
% flaglook: if 1 the signals are displayed.

%===============================================================
% the 3 following lines have to be removed if it is
% used as a function
directorydatafromIDC = '../../../AAdataI26calib/';
ihc      = 8;
flaglook = 1;

dirmat          = dir(sprintf('%ss%i/*.mat',directorydatafromIDC,ihc));
nbmats          = length(dirmat);
flagremainindex = zeros(nbmats,1);
Fs_Hz           = 20;
remainindex     = cell(nbmats,1);
for imatfile=63%1:nbmats,
    cdeload  = sprintf('load ''%ss%i/%s''',directorydatafromIDC,ihc,...
        dirmat(imatfile).name);
    eval(cdeload)
    Lrecords = length(records);
    N_ir     = zeros(Lrecords,1);
    maxcurr  = -inf;
    stdcurr  = -inf;
    for ir=1:Lrecords
        N_ir(ir) = length([records{ir}.data]);
        maxcurr  = max([maxcurr,max(abs([records{ir}.data]))]);
        stdcurr  = max([stdcurr,std([records{ir}.data])]);
    end
    if maxcurr<4*stdcurr
        flagremainindex(imatfile)=1;
    end
    
    if flaglook %and(length(NN)<=4,flaglook)
        setimesC_ihc  = zeros(nbmats,2);
        setimesH_ihc  = zeros(nbmats,2);
        signals       = zeros(34560000,2);
        idSc = 1;
        idSh = 1;
        cpC  = 1;
        cpH  = 1;
        
        for ir=1:Lrecords %min([Lrecords,4])
            if strcmp(records{ir}.channel,'BDF')
                switch records{ir}.station(4)
                    case 'C'
                        if cpC==1
                            cpC=cpC+1;
                            setimesC_ihc(imatfile,1) = records{ir}.stime;
                        end
                        auxC = records{ir}.etime;
                        LLC = length(records{ir}.data);
                        signalsC = [records{ir}.data];
                        signals(idSc:idSc+LLC-1,2)=signalsC;
                        idSc = idSc+LLC;
                    case 'H'
                        if cpH==1
                            cpH=cpH+1;
                            setimesH_ihc(imatfile,1) = records{ir}.stime;
                        end
                        auxH = records{ir}.etime;
                        LLH = length(records{ir}.data);
                        signalsH = [records{ir}.data];
                        signals(idSh:idSh+LLH-1,1)=signalsH;
                        idSh = idSh+LLH;
                end
            end
        end
        setimesC_ihc(imatfile,2) = auxC;
        setimesH_ihc(imatfile,2) = auxH;
        
        N = min([idSh,idSc]);
        signals = signals(1:N,:);
        subplot(211)
        plot((0:1000:N-1)/Fs_Hz/60,signals(1:1000:N,1) ,'r')
        subplot(212)
        plot((0:1000:N-1)/Fs_Hz/60,signals(1:1000:N,2) ,'b')
        drawnow
        pause
    end
    remainindex{imatfile} = find(flagremainindex==1);
end
