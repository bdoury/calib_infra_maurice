%====================== plotRsupdetails.m =============================
% Plot the details of the Rsups in the P frequency bands 
% as a function of times (along the 2 days).
% 
% The program uses the data extracted form the IDC and located in
% a given directory, here ../../../../AAsignals/
% and also the FB stored in 'filtercharacteristics/filtercharacteristics1.m'
%==========================================================================
clear
allcolors = ['b.';'r.';'m.';'c.';'g.';'k.';'rx';'yx';'mx';'rx';'kx';...
    'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];

addpath ../ZZtoolbox/

directorysignals  = '../../../../AAdataI26calib/';
% directory to save the figures
printdirectory    = ' ../../figures/';

%============== load the filter bank characteristics =======================
%  the useful variable is FILTERCHARACT
%======
% if generated by geneFB.m, use LOAD; if it is a .m program use RUN
% filtercharactfilename = 'filtercharacteristics/filtercharacteristics';
% cmdloadcharact        = sprintf('load(''%s'')',filtercharactfilename);
filtercharactfilename = '../filtercharacteristics/filtercharacteristics1.m';
cmdloadcharact        = sprintf('run(''%s'')',filtercharactfilename);
%======
eval(cmdloadcharact);
Pfilter = length(filtercharact);
%=============== Warning
% we have observed huge outliers in the following files:
% ihc==2, ifile== 33, i.e. sta2_Y2015_D219.mat from sample index 2.4e6
% ihc==8, ifile==63, i.e. sta8_Y2015_D280.mat from sample index 2.5e6
% ihc==5,ifile==62
% ihc==6,ifile==63
%??????????
%=====================
MSCthreshold = 0.98;
%=====================
ihc   = 1;%fix(5*rand)+1;
%===================== read data =========================
fileswithdotmat              = dir(sprintf('%ss%i/s%iyear*.mat',directorysignals,ihc,ihc));
nbmats                       = length(fileswithdotmat);
ifile                        = 31;%fix(nbmats*rand)+1;

allfrqsPfilters              = zeros(10000,1);
allRatioSupPfilters          = zeros(10000,1);
allSTDmodRatioSupPfilters    = zeros(10000,1);
allSTDphaseRatioSupPfilters  = zeros(10000,1);

allRatioInfPfilters          = zeros(10000,1);
allSTDmodRatioInfPfilters    = zeros(10000,1);
allSTDphaseRatioInfPfilters  = zeros(10000,1);
allmeanMSCcstPfilters        = zeros(10000,1);
nbofvaluesoverthreshold      = zeros(10000,1);


%==================================================
fullfilename_i      = fileswithdotmat(ifile).name;
dotlocation         = strfind(fullfilename_i,'.');
underscorelocation  = strfind(fullfilename_i,'_');
filenameonly        = fullfilename_i(setdiff(1:dotlocation-1,underscorelocation));
commandload         = sprintf('load %ss%i/%s',directorysignals,ihc,fullfilename_i);
eval(commandload)

if not(exist('directorysignals','var'))
    directorysignals    = '../AAsignals/';
end
Ts_sec = 1/Fs_Hz;

%============ Warning =======================
%============================================
if and(ihc==2,ifile==33)
    signals_centered = signals_centered(1:2.4e6,:);
end
if and(ihc==2,ifile==62)
    signals_centered = signals_centered([1:0.5e6 1e6:idScMin-1],:);
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

[SUTs, filteredsignals, allfrqsFFT_Hz, alltimes_sec, filterbank] = ...
    fbankanalysis(signals_centered,...
    filtercharact,Fs_Hz,MSCthreshold);
%============================================
P       = length(SUTs);
idipinf = zeros(P,1);
idipsup = zeros(P,1);
id1     = 1;
for ip=1:P
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
%%
constraintflag = 0;
displayhours = fix(size(signals_centered,1)/Fs_Hz/3600);
if constraintflag
    figure(1)
    clf
else
    figure(2)
    clf
end
    pos_ip = [0.07    0.1    0.9    0.1];
    deltay = 1.21*pos_ip(4);
    pos_time = pos_ip;
    pos_aux = get(gca,'position');
    pos_ip(2) = pos_aux(2)+deltay;

for ip=1:P
    [nbfrqs,nbtimeslots] = size(SUTs(ip).allMSCs.tabcst((idipinf(ip):idipsup(ip)),:));
    nbtimes = sum(sum(not(isnan((SUTs(ip).allMSCs.tabcst((idipinf(ip):idipsup(ip)),:))))));
    [nbtimes  nbfrqs*nbtimeslots]
    subplot('position',pos_ip);
    allfrqs_ip = allfrqsFFT_Hz{ip}(idipinf(ip):idipsup(ip));
    if constraintflag
        figure(1)
        if ip==1
            hpcolor=pcolor(alltimes_sec{ip}.SD/3600,allfrqs_ip(2:end),...
                20*log10(SUTs(ip).estimRsup.tabmodcst((idipinf(ip)+1:idipsup(ip)),:)));
                caxis([-20 1])

%             hpcolor=pcolor(alltimes_sec{ip}.SD/3600,allfrqs_ip,...
%                 (SUTs(ip).allMSCs.tabcst((idipinf(ip)+1:idipsup(ip)),:)));
        else
            hpcolor=pcolor(alltimes_sec{ip}.SD/3600,allfrqs_ip,...
                20*log10(SUTs(ip).estimRsup.tabmodcst((idipinf(ip):idipsup(ip)),:)));
                caxis([-20 1])

%              hpcolor=pcolor(alltimes_sec{ip}.SD/3600,allfrqs_ip,...
%                 (SUTs(ip).allMSCs.tabcst((idipinf(ip):idipsup(ip)),:)));           
        end
    else
        figure(2)
        if ip==1
            hpcolor=pcolor(alltimes_sec{ip}.SD/3600,allfrqs_ip(2:end),...
                20*log10(SUTs(ip).estimRsup.tabmod((idipinf(ip)+1:idipsup(ip)),:)));
        else

            hpcolor=pcolor(alltimes_sec{ip}.SD/3600,allfrqs_ip,...
                20*log10(SUTs(ip).estimRsup.tabmod((idipinf(ip):idipsup(ip)),:)));
        end
    end
    
%                  set(gca,'xlim',[10.5 11.4])

    set(gca,'xlim',[0 displayhours])
    caxis([-3 1])
    colorbar
    shading flat
    set(gca,'xtick',[])
    ylabel('Hz','fontname','times','fontsize',10)
    cmdqq=(sprintf('qq=[%3.2f;%3.2f];',allfrqs_ip([2 end-1])));
    eval(cmdqq)
    set(gca,'ytick',[qq(1),qq(2)])
    %     set(gca,'yscale','log')
    set(gca,'box','on')
    grid off
    pos_aux = get(gca,'position');
    pos_ip(2) = pos_aux(2)+deltay;
    set(gca,'fontname','times','fontsize',10)
    if ip==P
        if constraintflag
            title(sprintf('sensor: %i, year: %s, day:%s, threshold = %4.2f',...
                ihc,filenameonly(6:9),filenameonly(11:13),MSCthreshold),...
            'fontname','times','fontsize',12)
        else
        title(sprintf('sensor: %i, year: %s, day:%s',ihc,filenameonly(6:9),filenameonly(11:13)),...
            'fontname','times','fontsize',14)
        end
    end
    if constraintflag
        text(5,mean(allfrqs_ip([2 end-1])),sprintf('# events = %i',nbtimes))
    end
end
%===== time signals
pos_time(3) = pos_aux(3);
subplot('position',pos_time);
% maxsignals = max(max(abs(signals_centered(1:end,:))));
plot(alltimes_sec{ip}.signals(1:end)/3600,signals_centered(1:end,:))

%              set(gca,'xlim',[10.5 11.4])

             set(gca,'xlim',[0 displayhours])
set(gca,'ylim',[-12 12])
set(gca,'ytick',[])
set(gca,'fontname','times','fontsize',10)
xlabel('hours','fontname','times','fontsize',10)


HorizontalSize = 16;
VerticalSize   = 14;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');

if constraintflag
    fileprint = sprintf('%s2daysonIS26SUT%iyear%sday%saboveTH.eps',...
        printdirectory,ihc,filenameonly(6:9),filenameonly(11:13));
    figure(1)
else
    fileprint = sprintf('%s2daysonIS26SUT%iyear%sday%s.eps',...
        printdirectory,ihc,filenameonly(6:9),filenameonly(11:13));
    figure(2)
end

fileprintepscmd = sprintf('print -depsc -loose %s',fileprint);
fileeps2pdfcmd  = sprintf('!epstopdf %s',fileprint);
filermcmd       = sprintf('!rm %s',fileprint);

%
%         eval(fileprintepscmd)
%         eval(fileeps2pdfcmd)
%         eval(filermcmd)

%============================ END =========================================