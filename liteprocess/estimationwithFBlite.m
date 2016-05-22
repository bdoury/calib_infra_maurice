%====================== estimationwithFBlite.m ===========================
% Program estimates SUT response from the signals
% located in the directory directorysignals.
% The signals correspond to the pair of sensors
% SUT/SREF during a given duration T, typically 48 hours.
% Here we concatene NBCONCAT randomly chosen files.
%
% The evaluated parameters consist of the ratios, the STDs
% They are obtained by averaging on the period T.
% Results are plotted in figure 1
%=============
% Here we use the structure FILTERCHARACT and the list of frequencies.
% The only useful processing is the call to the function 
%              ESTIMSUTLITE.M
%=============  IMPORTANT: list of inputs to the developper ==============
% - filter.num and .den
% - allfreqsinfilter_Hz
% - Fs_Hz
% - the signals (signals_centered variable in the following)
% - ref sensor response
%=========================================================================
clear

% the following lines can be changed by the user:
MSCthreshold     = 0.9;
FLAGsavesmall    = 0;
Fs_Hz            = 20;
ihc              = 1;
trimpercent      = 1;
nbrandomconcat   = 1;
frequencylist_Hz = logspace(-2,log10(6),30);
%===
% directories
directorysignals    = '../../../AAdataI26calib/';
directoryFB         = '../fullprocess/filtercharacteristics/';
%=========================================================================
%=========================================================================
%====== load the filter bank characteristics
%  the useful structure is FILTERCHARACT
filtercharactfilename = 'filtercharacteristics1.m';
cmdloadfilter         = sprintf('run(''%s%s'')',directoryFB,...
    filtercharactfilename);
eval(cmdloadfilter);
%===================== read data =========================
fileswithdotmat        = dir(sprintf('%ss%i/s%iy*.mat',...
    directorysignals,ihc,ihc));
nbmats                 = length(fileswithdotmat);
signals                = [];
indperm                = 1:nbmats; %randperm(nbmats);
alldates               = cell(nbrandomconcat,1);
selectedlist           = (1:nbrandomconcat);%   indperm(1:nbrandomconcat);
for indfile = 1:nbrandomconcat
    ifile              = selectedlist(indfile);
    fullfilename_i     = fileswithdotmat(ifile).name;
    dotlocation        = strfind(fullfilename_i,'.');
    underscorelocation = strfind(fullfilename_i,'_');
    filenameonly       = fullfilename_i(...
        setdiff(1:dotlocation-1,...
        underscorelocation));
    commandload        = sprintf('load %ss%i/%s',...
        directorysignals,ihc,fullfilename_i);
    eval(commandload)
    aux = str2double(fullfilename_i(21:22));
    if aux<9
        straux = ['0' num2str(aux+1)];
    else
        straux = num2str(aux+1);
    end
    date_i             = ...
        sprintf('%s/%s/%s-%s',fullfilename_i(7:10),...
        fullfilename_i(16:17),fullfilename_i(21:22),...
        straux);
    alldates{indfile} = date_i;
    signals           = [signals;signals_centered];
end
disp('************************************************')
sortalldates = sort(alldates);
display(sprintf('Station %i:',ihc));
display(sort(alldates))
%%
txtlatex = cell(nbrandomconcat,1);
for ii=1:nbrandomconcat
   aux = sprintf('%stt %s %s','{\',[sortalldates{ii}],'}\\');
   txtlatex{ii} = aux;
end
cell2mat(txtlatex)
%%
disp('************* start process ********************')
%===============================================================
%===============================================================
%=============== processing function call ======================
%===============================================================
tic
[Rsup, freqslin, STDmoduleR, STDphaseR, nboverTH, R] = ...
    estimSUTlite ...
    (signals, filtercharact, frequencylist_Hz, ...
    Fs_Hz, MSCthreshold, trimpercent);
toc
%%
%===============================================================
%============================ for plotting =====================
%===============================================================
%===============================================================
figure(1)
semilogx(freqslin, 20*log10(abs(Rsup)),...
    'or','markerfacecolor','r')
hold on 
semilogx(frequencylist_Hz, 1.5*ones(length(frequencylist_Hz),1),...
    'ob','markerfacecolor','b')
hold off
set(gca,'fontname','times','fontsize',12)
grid on
set(gca,'xlim',[0.01 7],'ylim',[-1 2])
hold off
title(sprintf('%i days on station %i',2*nbrandomconcat,ihc),...
    'fontname','times','fontsize',14)
xlabel('frequency - Hz')
ylabel('Rsup - dB')
%==
HorizontalSize = 12;
VerticalSize   = 10;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');

%===============================================================
%============================ for printing =====================
%===============================================================
numexample = 3;
printdirectory  = ' ../calibtextelite/';
fileprint = sprintf('%sexample%ionstation%i.eps',...
    printdirectory,numexample,ihc);
figure(1)
fileprintepscmd = sprintf('print -depsc -loose %s',fileprint);
fileeps2pdfcmd  = sprintf('!epstopdf %s',fileprint);
filermcmd       = sprintf('!rm %s',fileprint);


%     eval(fileprintepscmd)
%     eval(fileeps2pdfcmd)
%     eval(filermcmd)