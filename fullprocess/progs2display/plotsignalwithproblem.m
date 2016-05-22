%====================== plotsignalwithproblem.m ==========================
% This program plots the signals and some problems
% inputs are the directory of analyzed signals,
% the station number and the date 
%===========================================

clear
addpath ../ZZtoolbox/
%== signal directory
directorysignals  = '../../../../AAdataI26calib/';
%== station number
ihc               = 1;
%== selected date
date_select       = '2015/07/09';

%============ examples ======================
%=========== non exhaustive list ============
% if and(ihc==1,ifile==62)
% if and(ihc==2,ifile==33)
%     signals_centered = signals_centered(1:2.4e6,:);
% end
% if and(ihc==2,ifile==67)
%     signals_centered = signals_centered([1:0.5e6 1e6:idScMin-1],:);
% end
% if and(ihc==5,ifile==62)
%     signals_centered = signals_centered(1:2.4e6,:);
% end
% if and(ihc==6,ifile==63)
%     signals_centered = signals_centered(1:2.14e6,:);
% end
% if and(ihc==8,ifile==63)
%     signals_centered = signals_centered(1:2.5e6,:);
% end
%=========================================================================
fileswithdotmat         = dir(sprintf('%ss%i/s%iy*.mat',...
    directorysignals,ihc,ihc));
nbmats                  = length(fileswithdotmat);
flagnofile_i            = 1;
for ifile=1:nbmats
    fullfilename_i      = fileswithdotmat(ifile).name;
    dotlocation         = strfind(fullfilename_i,'.');
    underscorelocation  = strfind(fullfilename_i,'_');
    filenameonly        = fullfilename_i(...
        setdiff(1:dotlocation-1,underscorelocation));
    date_i = sprintf('%s/%s/%s',fullfilename_i(7:10),...
        fullfilename_i(16:17),fullfilename_i(21:22));
    if strcmp(date_i,date_select)
        flagnofile_i = 0;
        ifile_select = ifile;
        fullfilename_select = fullfilename_i;
    end
end
fullfilename_i      = fileswithdotmat(ifile_select).name;
dotlocation         = strfind(fullfilename_i,'.');
underscorelocation  = strfind(fullfilename_i,'_');
filenameonly        = fullfilename_i(setdiff(1:dotlocation-1,...
    underscorelocation));
commandload         = sprintf('load %ss%i/%s',directorysignals,...
    ihc,fullfilename_i);
date_ii = sprintf('%s/%s/%s',filenameonly(7:10),...
    filenameonly(16:17),filenameonly(21:22));
eval(commandload)
subplot(311)
plot(signals_centered(:,1))
subplot(312)
plot(signals_centered(:,2))
subplot(313)
plot(signals_centered)
for is=1:2
    subplot(2,3,3*is-2)
    plot((0:size(signals_centered,1)-1)/Fs_Hz/3600,...
        signals_centered(:,is))
    set(gca,'xlim',[0 size(signals_centered,1)/Fs_Hz/3600])
    set(gca,'fontname','times','fontsize',10)
    xlabel('hours','fontname','times','fontsize',10)
    if is==1
        title(sprintf('%s on H1',date_ii),'fontname','times','fontsize',8)
    else
        title(sprintf('%s on C1',date_ii),'fontname','times','fontsize',8)
    end
end
for is=1:2
    subplot(2,3,3*is-1)
    id1=8640*60;id2=8880*60;
    plot((0:id2-id1)/Fs_Hz/60,signals_centered(id1:id2,is))
    set(gca,'xlim',[0 (id2-id1)/Fs_Hz/60])
    %             set(gca,'xticklabel',[])
    if is ==1
        set(gca,'fontname','times','fontsize',10)
        xlabel('minutes','fontname','times','fontsize',10)
        title('Zoom on the first burst','fontname','times','fontsize',8)
    else
        set(gca,'fontname','times','fontsize',10)
        xlabel('minutes','fontname','times','fontsize',10)
        title('Zoom on the second burst','fontname','times','fontsize',8)
    end
end
for is=1:2
    subplot(2,3,3*is)
    plot((0:id2-id1)/Fs_Hz/60,signals_centered((id1:id2)+100000,is))
    set(gca,'xlim',[0 (id2-id1)/Fs_Hz/60])
    set(gca,'xticklabel',[])
    set(gca,'ylim',[-0.5 0.5])
    set(gca,'fontname','times','fontsize',10)
    xlabel('minutes','fontname','times','fontsize',10)
    title('Zoom outside of the bursts','fontname','times','fontsize',8)
end

set(gca,'fontname','times','fontsize',10)

HorizontalSize = 14;
VerticalSize   = 12;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
%     set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');

printdirectory  = ' ../../figures/';

date_ii = sprintf('%s/%s/%s',filenameonly(7:10),...
    filenameonly(16:17),filenameonly(21:22));


fileprint = sprintf('%sclickdetails%iyear%smonth%sday%saboveTH.eps',...
    printdirectory,ihc,filenameonly(7:10),filenameonly(16:17),...
    filenameonly(21:22));

fileprintepscmd = sprintf('print -depsc -loose %s',fileprint);
fileeps2pdfcmd  = sprintf('!epstopdf %s',fileprint);
filermcmd       = sprintf('!rm %s',fileprint);

%
%     eval(fileprintepscmd)
%     eval(fileeps2pdfcmd)
%     eval(filermcmd)
%============================ END =======================================
