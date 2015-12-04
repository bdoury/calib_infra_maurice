%====================== estimationwithFB.m ===============================
% Program evaluates some parameters from the signals
% located in the directory directorysignals. 
% The signals correspond to the pair of sensors 
% SUT/SREF during a given duration T, typically 48 hours.
% 
% The evaluated parameters consist of the ratios, the STDs, etc
% They are obtained by averaging on the period T.
% They are saved in directoryresults.
% 
% Then for display use programs in the directory
% progs2display.
%
%=========================================================================
% used function
%     - fbankanalysis.m
%
%=========================================================================
% Rk: 
% The 3 following quantities do not depend on the
% index ifile. Therefore they could be performed
% outside of the loop on ifile. Indeed
% SUTs(ip).indexinsidefreqband can be performed outside
% of the function fbankanalysis.m
% idipinf       = zeros(Pfilter,1);
% idipsup       = zeros(Pfilter,1);
% cumsumnbfq_ip = zeros(Pfilter,2);
% In the following program they are performed at each step of the loop
% using a quatity performed by the function fbankanalysis.m
%=========================================================================
%
%==================== Warning ===============
% we have observed huge outliers/problems in the following files:
% ihc==1, date==2015/08/09 , shift between the 2 signals for the full
% duration
% ihc==1, date==2015/10/07 ,from sample index 2.4e6
% ihc==1, date==2015/10/09 ,from sample index 2.4e6
% ihc==2, date==2015/08/07 ,from sample index 2.4e6
% ihc==2, date==2015/10/05 ,from sample index 2.4e6
% ihc==5, date==2015/10/05
% ihc==6, date==2015/10/07
% ihc==8, date==2015/10/07 , from sample index 2.5e6
% in this version we remove them
%
%=========================================================================

clear
addpath ZZtoolbox/

%=====================
MSCthreshold     = 0.98;
FLAGsaveall      = 0;
FLAGsavesmall    = 0;
nofilterbankflag = 0;
flagtheoSTD      = 0;

%=====================
%=== directory of input signals
directorysignals    = '../../../AAdataI26calib/';
%=== directory of output results (huge records)
% if FLAGsaveall=1
directoryresultsALL = 'BBresults'; 
% if FLAGsavesmall=1
directoryresults    = sprintf('AAresultswithFB%i_10',fix(MSCthreshold*100));

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

if nofilterbankflag
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

for ihc = 1 %:8,  ihc is the site - index of H and C
    %===================== read data =========================
    fileswithdotmat              = dir(sprintf('%ss%i/s%iy*.mat',...
        directorysignals,ihc,ihc));
    nbmats                       = length(fileswithdotmat);
    nbmats=1;   % select nbmats period of 2 days out of the whole set
    ifileselect = 20;   % index (which nbmats period of 2 days)
    
    %====== Useful evaluated parameters for general purposes
    % 10000 means that we dont  know here the values. It is 
    % adjusted at the end of the loop.
    allfrqsPfilters              = zeros(10000,nbmats);
    allRatioSupPfilters          = zeros(10000,nbmats);
    allSTDmodRatioSupPfilters    = zeros(10000,nbmats);
    allSTDphaseRatioSupPfilters  = zeros(10000,nbmats);
    theoreticalSTDmod            = zeros(10000,nbmats);
    theoreticalSTDphase_rad      = zeros(10000,nbmats);
    
    % the RatioInf associated parameters are not considered
    % in the final use.
    allRatioInfPfilters          = zeros(10000,nbmats);
    allSTDmodRatioInfPfilters    = zeros(10000,nbmats);
    allSTDphaseRatioInfPfilters  = zeros(10000,nbmats);
    allmeanMSCcstPfilters        = zeros(10000,nbmats);
    nbofvaluesoverthreshold      = zeros(10000,nbmats);
    allScpPfilters               = zeros(3,10000,nbmats);
    %==================================================
    for ifile=1:nbmats, ifile,
        tic
        fullfilename_i      = fileswithdotmat(ifileselect).name;
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

        date_i = sprintf('%s/%s/%s',fullfilename_i(7:10),...
            fullfilename_i(16:17),fullfilename_i(21:22));

        %============ Warning =======================
        %============================================
        %============================================
        %== some manual checks 
        if and(ihc==1,date_i=='2015/10/07')
            signals_centered = signals_centered(0.6e6:2.2e6,:);
        end
        if and(ihc==1,date_i=='2015/10/09')
            signals_centered = signals_centered(0.6e6:Ntotal,:);
        end
        if and(ihc==2,date_i=='2015/08/07')
            signals_centered = signals_centered(1:2.4e6,:);
        end
        if and(ihc==2,date_i=='2015/10/05')
            signals_centered = ...
                signals_centered([1:0.5e6 1e6:Ntotal],:);
        end
        if and(ihc==2,date_i=='2015/10/13')
            signals_centered = ...
                signals_centered([1e6:2.2e6 2.3e6:Ntotal],:);
        end
        
        if and(ihc==5,date_i=='2015/10/05')
            signals_centered = signals_centered(1:2.4e6,:);
        end
        
        if and(ihc==6,date_i=='2015/10/07')
            signals_centered = signals_centered(1:2.14e6,:);
        end
        
        if and(ihc==8,date_i=='2015/10/07')
            signals_centered = signals_centered(1:2.5e6,:);
        end
        %============================================        
        %============================================
        % notice that the SUTs is not saved, therefore we have only the
        % last associated the laxt index which is NBMATS
        
        [SUTs, filteredsignals, allfrqsFFT_Hz, ...
            alltimes_sec, filterbank] = ...
            fbankanalysis(signals_centered,...
            filtercharact,Fs_Hz,MSCthreshold,flagtheoSTD);
        %============================================
        % These three quantities, idipinf, idipsup and cumsumnbfq_ip
        % do not depend on the index ifile and do depend only on the
        % filtercharac parameters. Therefore they could be performed
        % outside of the loop on ifile. Indeed the variable 
        % SUTs(ip).indexinsidefreqband, provided by the function
        % fbankanalysis can be performed outside
        % 
        idipinf       = zeros(Pfilter,1);
        idipsup       = zeros(Pfilter,1);
        cumsumnbfq_ip = zeros(Pfilter,2);
        id1     = 1;
        for ip=1:Pfilter
            cumsumnbfq_ip(ip,1) = id1;
            idipinf(ip) = SUTs(ip).indexinsidefreqband(1);
            idipsup(ip) = SUTs(ip).indexinsidefreqband(2);
            id2         = id1+(idipsup(ip)-idipinf(ip));
            cumsumnbfq_ip(ip,2) = id2;
            id1         = id2+1;
        end
        %====== Useful evaluated parameters for general purposes
        for ip=1:Pfilter
            id1               = cumsumnbfq_ip(ip,1);
            id2               = cumsumnbfq_ip(ip,2);
            allRatioSupPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRsup.modcst(idipinf(ip):idipsup(ip)) .* ...
                exp(1i*SUTs(ip).estimRsup.phasecst(idipinf(ip):idipsup(ip)));
            allSTDmodRatioSupPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRsup.stdmodcst(idipinf(ip):idipsup(ip));
            allSTDphaseRatioSupPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRsup.stdphasecst(idipinf(ip):idipsup(ip));
            
            allRatioInfPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRinf.modcst(idipinf(ip):idipsup(ip)) .* ...
                exp(1i*SUTs(ip).estimRinf.phasecst(idipinf(ip):idipsup(ip)));
            allSTDmodRatioInfPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRinf.stdmodcst(idipinf(ip):idipsup(ip));
            allSTDphaseRatioInfPfilters(id1:id2,ifile) = ...
                SUTs(ip).estimRinf.stdphasecst(idipinf(ip):idipsup(ip));
            
            allmeanMSCcstPfilters(id1:id2,ifile) = ...
                nanmean(SUTs(ip).allMSCs.tabcst(idipinf(ip):idipsup(ip),:),2);
            allfrqsPfilters(id1:id2,ifile) = ...
                SUTs(ip).frqsFFT_Hz(idipinf(ip):idipsup(ip))';
            nbofvaluesoverthreshold(id1:id2,ifile) = ...
                sum(not(isnan(SUTs(ip).allMSCs.tabcst(idipinf(ip):idipsup(ip),:))),2);
            allScpPfilters(:,id1:id2,ifile) = ...
                squeeze(SUTs(ip).spectralmatrix(:,idipinf(ip):idipsup(ip)));
            
            theoreticalSTDmod(id1:id2,ifile) = SUTs(ip).theoSTDmodonRsup;
            theoreticalSTDphase_rad(id1:id2,ifile) = ...
                SUTs(ip).theoSTDphase_rd;
        end
%         if FLAGsaveall
%             comsave = ...
%                 sprintf('save %s/s%i/resultssta26sensor%s', ...
%                 directoryresultsALL,ihc,filenameonly);
%             eval(comsave);
%         end
         toc
    end
    
    allRatioSupPfilters         = allRatioSupPfilters(1:id2,:);
    allSTDmodRatioSupPfilters   = allSTDmodRatioSupPfilters(1:id2,:);
    allSTDphaseRatioSupPfilters = allSTDphaseRatioSupPfilters(1:id2,:);
    
    allRatioInfPfilters         = allRatioInfPfilters(1:id2,:);
    allSTDmodRatioInfPfilters   = allSTDmodRatioInfPfilters(1:id2,:);
    allSTDphaseRatioInfPfilters = allSTDphaseRatioInfPfilters(1:id2,:);
    
    allfrqsPfilters             = allfrqsPfilters(1:id2,1);
    allmeanMSCcstPfilters       = allmeanMSCcstPfilters(1:id2,:);
    nbofvaluesoverthreshold     = nbofvaluesoverthreshold(1:id2,:);
    
    allScpPfilters              = allScpPfilters(:,1:id2,:);
    if flagtheoSTD
        theoreticalSTDmod           = theoreticalSTDmod(1:id2,:);
        theoreticalSTDphase_rad     = theoreticalSTDphase_rad(1:id2,:);
    end
    if FLAGsavesmall
        comsave = ...
            sprintf('save %s/resultssta26sensor%i',...
            directoryresults,ihc);
        clear signals_centered
        clear filteredsignals
        clear alltimes_sec
        clear SUTs
        eval(comsave);
    end
end
%============================ END =========================================