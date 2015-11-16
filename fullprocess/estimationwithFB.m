%====================== estimationwithFB.m ===============================
% program evaluates some parameters from the signals
% located in directorysignals. These parameters are 
% saved in directoryresults.
% Then for display use programs in the directory
% progs2display.
%
%=========================================================================
% used function
%     - fbankanalysis.m
%
%=========================================================================

%=============== Warning ===============
% we have observed huge outliers in the following files:
% ihc==2, ifile== 33, i.e. sta2_Y2015_D219.mat from sample index 2.4e6
% ihc==8, ifile==63, i.e. sta8_Y2015_D280.mat from sample index 2.5e6
% ihc==5,ifile==62
% ihc==6,ifile==63
% in this version we remove them but in some saved data that has not be
% done
%=========================================================================

clear
allcolors = ['b.';'r.';'m.';'c.';'g.';'k.';'rx';'yx';'mx';'rx';'kx';...
    'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];
%=====================
MSCthreshold = 0.98;
%=====================

FLAGsaveall   = 0;
FLAGsavesmall = 1;
addpath ZZtoolbox/

%=== directory of input signals
directorysignals      = '../../../AAdataI26calib/';
%=== directory of output results
directoryresultsALL   = 'BBresults'; % if FLAGsaveall=1
directoryresults      = sprintf('AAresultswithFB%i',fix(MSCthreshold*100));% if FLAGsavesmall=1

%============== load the filter bank characteristics =====================
%  the useful variable is FILTERCHARACT
%======
% if generated by geneFB.m, use LOAD; if it is a .m program use RUN
% filtercharactfilename = 'filtercharacteristics/filtercharacteristics';
%cmdloadfilter        = sprintf('load(''%s'')',filtercharactfilename);
filtercharactfilename = 'filtercharacteristics/filtercharacteristics1.m';
cmdloadfilter         = sprintf('run(''%s'')',filtercharactfilename);
eval(cmdloadfilter);

%======
nofilterflag = 0;

if nofilterflag
    clear filtercharact;
    filtercharact(1).designname     = 'butter';
    filtercharact(1).Norder         = 2;
    filtercharact(1).Wlow_Hz        = 0.01;
    filtercharact(1).Whigh_Hz       = 5;
    filtercharact(1).SCPperiod_sec  = 300;
    filtercharact(1).windowshape    = 'hann';
    filtercharact(1).overlapDFT     = 0.5;
    filtercharact(1).overlapSCP     = 0;
    filtercharact(1).ratioDFT2SCP   = 5;
end

%=====================
Pfilter = length(filtercharact);

% if and(Pfilter==1, filtercharact(Pfilter).Norder==0)
%     filtercharact(Pfilter).Wlow_Hz  = 0.001;
%     filtercharact(Pfilter).Whigh_Hz = 10;
% end

for ihc = 6:8, ihc
    %=====================
    % under test = 1, reference = 2
    %===================== read data =========================
    fileswithdotmat              = dir(sprintf('%ss%i/s%iy*.mat',...
        directorysignals,ihc,ihc));
    nbmats                       = length(fileswithdotmat);

    allfrqsPfilters              = zeros(10000,nbmats);
    allRatioSupPfilters          = zeros(10000,nbmats);
    allSTDmodRatioSupPfilters    = zeros(10000,nbmats);
    allSTDphaseRatioSupPfilters  = zeros(10000,nbmats);
    
    allRatioInfPfilters          = zeros(10000,nbmats);
    allSTDmodRatioInfPfilters    = zeros(10000,nbmats);
    allSTDphaseRatioInfPfilters  = zeros(10000,nbmats);
    allmeanMSCcstPfilters        = zeros(10000,nbmats);
    nbofvaluesoverthreshold      = zeros(10000,nbmats);
%     RsupmodtabNoThreshold        = cell(nbmats,Pfilter);
%     MSCtabNoThreshold            = cell(nbmats,Pfilter);
    spectralmatrixtab            = cell(nbmats,Pfilter);
    %==================================================
    for ifile=1:nbmats, ifile,tic
        fullfilename_i      = fileswithdotmat(ifile).name;
        dotlocation         = strfind(fullfilename_i,'.');
        underscorelocation  = strfind(fullfilename_i,'_');
        filenameonly        = fullfilename_i(...
            setdiff(1:dotlocation-1,...
            underscorelocation));
        commandload         = sprintf('load %ss%i/%s',...
            directorysignals,ihc,fullfilename_i);
        eval(commandload)
        Ts_sec = 1/Fs_Hz;
        Ntotal = size(signals_centered,1);

        %============ Warning =======================
        %============================================
        %============================================
        %== some manual checks 
        if and(ihc==2,ifile==33)
            signals_centered = signals_centered(1:2.4e6,:);
        end
        
        if and(ihc==2,ifile==62)
            signals_centered = ...
                signals_centered([1:0.5e6 1e6:Ntotal],:);
        end
        
        if and(ihc==5,ifile==62)
            signals_centered = signals_centered(1:2.4e6,:);
        end
        
        if and(ihc==6,ifile==63)
            signals_centered = signals_centered(1:2.14e6,:);
        end
        
        if and(ihc==8,ifile==63)
            signals_centered = signals_centered(1:2.5e6,:);
        end
        %============================================        
        %============================================
        % notice that the SUTs is not saved, therefore we have only the
        % last associated the laxt index which is NBMATS
        
        [SUTs, filteredsignals, allfrqsFFT_Hz, ...
            alltimes_sec, filterbank] = ...
            fbankanalysis(signals_centered,...
            filtercharact,Fs_Hz,MSCthreshold);
        %============================================
        idipinf = zeros(Pfilter,1);
        idipsup = zeros(Pfilter,1);
        id1     = 1;
        for ip=1:Pfilter
            idipinf(ip) = SUTs(ip).indexinsidefreqband(1);
            idipsup(ip) = SUTs(ip).indexinsidefreqband(2);
            id2         = id1+(idipsup(ip)-idipinf(ip));
            allRatioSupPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRsup.modcst(idipinf(ip):idipsup(ip)) .* ...
                exp(1i*SUTs(ip).estimRsup.phasecst(idipinf(ip):idipsup(ip)));
            allSTDmodRatioSupPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRsup.stdmodcst(idipinf(ip):idipsup(ip));
            allSTDphaseRatioSupPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRsup.phasecst(idipinf(ip):idipsup(ip));
            
            allRatioInfPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRinf.modcst(idipinf(ip):idipsup(ip)) .* ...
                exp(1i*SUTs(ip).estimRinf.phasecst(idipinf(ip):idipsup(ip)));
            allSTDmodRatioInfPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRinf.stdmodcst(idipinf(ip):idipsup(ip));
            allSTDphaseRatioInfPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRinf.phasecst(idipinf(ip):idipsup(ip));
            
            allmeanMSCcstPfilters(id1:id2,ifile) = ...
                nanmean(SUTs(ip).allMSCs.tabcst(idipinf(ip):idipsup(ip),:),2);
            allfrqsPfilters(id1:id2,ifile) = ...
                SUTs(ip).frqsFFT_Hz(idipinf(ip):idipsup(ip))';
            nbofvaluesoverthreshold(id1:id2,ifile) = ...
                sum(not(isnan(SUTs(ip).allMSCs.tabcst(idipinf(ip):idipsup(ip),:))),2);
            id1 = id2+1;
            
%             RsupmodtabNoThreshold{ifile,ip} = ...
%                 squeeze(SUTs(ip).estimRsup.tabmod(idipinf(ip):idipsup(ip),:));
%             MSCtabNoThreshold{ifile,ip} = ...
%                 squeeze(SUTs(ip).allMSCs.tab(idipinf(ip):idipsup(ip),:));
            spectralmatrixtab{ifile,ip} = ...
                squeeze(SUTs(ip).spectralmatrix(:,idipinf(ip):idipsup(ip)));
        end
%         if FLAGsaveall
%             comsave = ...
%                 sprintf('save %s/s%i/resultssta26sensor%s', ...
%                 directoryresultsALL,ihc,filenameonly);
%             eval(comsave);
%         end
        toc
    end
    
    allRatioSupPfilters         = allRatioSupPfilters(1:id1-1,:);
    allSTDmodRatioSupPfilters   = allSTDmodRatioSupPfilters(1:id1-1,:);
    allSTDphaseRatioSupPfilters = allSTDphaseRatioSupPfilters(1:id1-1,:);
    
    allRatioInfPfilters         = allRatioInfPfilters(1:id1-1,:);
    allSTDmodRatioInfPfilters   = allSTDmodRatioInfPfilters(1:id1-1,:);
    allSTDphaseRatioInfPfilters = allSTDphaseRatioInfPfilters(1:id1-1,:);
    
    allfrqsPfilters             = allfrqsPfilters(1:id1-1,1);
    allmeanMSCcstPfilters       = allmeanMSCcstPfilters(1:id1-1,:);
    nbofvaluesoverthreshold     = nbofvaluesoverthreshold(1:id1-1,:);

    if FLAGsavesmall
        comsave = ...
            sprintf('save %s/resultssta26sensor%i',...
            directoryresults,ihc+100);
        clear signals_centered
        clear filteredsignals
        clear alltimes_sec
        clear SUTs
        eval(comsave);
    end
end
%============================ END =========================================