%====================== plotsignalclick.m =============================
clear
allcolors = ['b.';'r.';'m.';'c.';'g.';'k.';'rx';'yx';'mx';'rx';'kx';...
    'c.';'k.';'r.';'c.';'m.';'g.';'b.';'k.';'r.';'c.';'m.';'g.';'k.'];

addpath ../ZZtoolbox/
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
%============ Warning =======================
%============================================
%============================================
%============================================
%============================================
% if and(ihc==1,ifile==63)
% if and(ihc==2,ifile==33)
%     signals_centered = signals_centered(1:2.4e6,:);
% end
% if and(ihc==2,ifile==62)
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
%============================================
%============================================
%============================================
%=====================
directorysave2daysignals  = '../../../../AAdataI26calib/';
%===================== read data =========================
% nbmats                       = length(fileswithdotmat);

%==================================================
for ihc   = 2
    fileswithdotmat     = dir(sprintf('%ss%i/s%iy*.mat', ...
        directorysave2daysignals,ihc,ihc));
    nbmats = length(fileswithdotmat);
    for ifile = 67%1:nbmats,ifile
        
        fullfilename_i      = fileswithdotmat(ifile).name;
        dotlocation         = strfind(fullfilename_i,'.');
        underscorelocation  = strfind(fullfilename_i,'_');
        filenameonly        = fullfilename_i(setdiff(1:dotlocation-1,underscorelocation));
        commandload         = sprintf('load %ss%i/%s',directorysave2daysignals,...
            ihc,fullfilename_i);
        eval(commandload)
        for is=1:2
            subplot(2,3,3*is-2)
            plot((0:size(signals_centered,1)-1)/Fs_Hz/60,signals_centered(:,is))
            set(gca,'xlim',[0 size(signals_centered,1)/Fs_Hz/60])
                set(gca,'fontname','times','fontsize',10)
                xlabel('hours','fontname','times','fontsize',10)
            title(fullfilename_i,'fontname','times','fontsize',8)

        end
        for is=1:2
            subplot(2,3,3*is-1)
            id1=8640*60;id2=8880*60;
            plot((0:id2-id1)/Fs_Hz/60,signals_centered(id1:id2,is))
            set(gca,'xlim',[0 (id2-id1)/Fs_Hz/60])
%             set(gca,'xticklabel',[])
                set(gca,'fontname','times','fontsize',10)
                xlabel('minutes','fontname','times','fontsize',10)
                    title('Zoom on the first burst','fontname','times','fontsize',8)
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
    end
    


    set(gca,'fontname','times','fontsize',10)
    
    
    HorizontalSize = 12;
    VerticalSize   = 8;
    set(gcf,'units','centimeters');
    set(gcf,'paperunits','centimeters');
    set(gcf,'PaperType','a3');
%     set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
    set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
    set(gcf,'color', [1,1,0.92]);
    set(gcf, 'InvertHardCopy', 'off');
    
    printdirectory  = ' ../../figures/';
    
    fileprint = sprintf('%sclickdetails%iyear%sday%saboveTH.eps',...
        printdirectory,ihc,filenameonly(7:10),filenameonly(14:16));
    
    fileprintepscmd = sprintf('print -depsc -loose %s',fileprint);
    fileeps2pdfcmd  = sprintf('!epstopdf %s',fileprint);
    filermcmd       = sprintf('!rm %s',fileprint);
    
    %
%             eval(fileprintepscmd)
%             eval(fileeps2pdfcmd)
% `            eval(filermcmd)
end
%============================ END =========================================


%%
figure(2)

subplot(211)
plot((0:size(signals_centered,1)-1)/Fs_Hz/60,signals_centered(:,1))

subplot(212)
plot((0:size(signals_centered,1)-1)/Fs_Hz/60,signals_centered(:,2))

