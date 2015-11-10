function remainindex = lookatdata(directorydatafromIDC,ihc, flaglook)

clear
directorydatafromIDC  = '../../../AAdataI26calib/';
ihc = 2;
flaglook = 0;

dirmat = dir(sprintf('%ss%i/*.mat',directorydatafromIDC,ihc));
nbmats = length(dirmat);
flagremainindex = zeros(nbmats,1);
Fs_Hz = 20;


for imatfile=1:nbmats,
    cdeload  = sprintf('load ''%ss%i/%s''',directorydatafromIDC,ihc,...
        dirmat(imatfile).name);
    eval(cdeload)
    Lrecords = length(records);
    NN       = zeros(Lrecords,1);
    maxcurr  = -inf;
    for ir=1:Lrecords
        NN(ir) = length([records{ir}.data]);
        maxcurr = max([maxcurr,max(abs([records{ir}.data]))]);
    end
    if and(length(NN)<=4, maxcurr<300)
        flagremainindex(imatfile)=1;
    end
    
    if and(length(NN)<=4,flaglook)
        setimesC_ihc  = zeros(nbmats,2);
        setimesH_ihc  = zeros(nbmats,2);
        signals       = zeros(34560000,2);
        idSc = 1;
        idSh = 1;
        cpC  = 1;
        cpH  = 1;
        
        for ir=1:min([length(NN),4])
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
end
remainindex = find(flagremainindex==1);