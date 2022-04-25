function graph_loc
   %% Graph the lines of code
   datafile; dat = loc; t=dat(:,5);
   clf
   hold on;
   hh(1)=plot(t,1e-3*dat(:,1),            '-','Color',[2,0,0]/3,'LineWidth',3);
   hh(2)=plot(t,1e-3*dat(:,1:3)*[1;-1;-1],'-','Color',[1,1,0]/3,'LineWidth',3);
   hh(3)=plot(t,1e-3*dat(:,2),            '-','Color',[0,2,0]/3,'LineWidth',3);
   hh(4)=plot(t,1e-3*dat(:,3),            '-','Color',[0,0,2]/3,'LineWidth',3);

   k=1;for yy=2004:2022;
     tt(k) = (datenum(yy,1,1)-datenum(1970,1,1))*86400;
     tl{k} = num2str(yy); if rem(yy,2)==0; tl{k}=''; end
   k=k+1; end
   set(gca,'xtick',tt)
   set(gca,'xticklabel',tl);
   ylim([0, ylim*[0;1]]);
   xlim([min(t),max(t)+0e7]);
   releaselines
%  uistack(hh,'top'); %%% Screws up colours
   hold off;
   set(gca,'Box','off')
   set(gcf,'PaperPosition',[1,1,6,2.5]);
   legend('Total','Eidors','Tutorials','Dev','location','northwest');
   ylabel('Lines of Code (x1000)')
   print -dpdf fig_loc.pdf
   !LD_LIBRARY_PATH="" pdfcrop fig_loc.pdf fig_loc.pdf

function releaselines

release= [
(datenum(2019,06,30)-datenum(1970,1,1))*86400; % eidors-v3.9.1
(datenum(2018,06,1)-datenum(1970,1,1))*86400; % eidors-v3.9.1
(datenum(2017,06,21)-datenum(1970,1,1))*86400; % eidors-v3.9
1431302400; %eidors-v3.8     2015-05-11
1369785600; %eidors-v3.7.1   2013-05-29
1366675200; %eidors_v3.7     2013-04-23
%1366243200; %eidors-v3.7rc   2013-04-18
1341446400; %eidors-v3.6     2012-07-05
1310688000; %eidors-v3.5     2011-07-15
1278633600; %eidors-v3.4     2010-07-09
1216771200; %eidors-v3.3     2008-07-23
%1213315200; %eidors-v3.3rc1  2008-06-13
1188518400; %EIDORS 3.2      2007-08-31
1156204800; %EIDORS 3.1      2006-08-22
1138060800; %EIDORS 3.1RC1   2006-01-24
1130803200; %EIDORS 3.0      2005-11-01
%1129680000; %EIDORS 3.0RC1   2005-10-19
];
   hh=plot([1;1]*release(:)',[0,ylim*[0;1]], 'Color',[2,2,2]/4,'LineWidth',2);
%  uistack(hh,'bottom'); % screws up legend
   set(hh,'HandleVisibility','off'); % Don't show up in legend
   hh= get(gca,'Children'); 
   set(gca,'Children',[hh(end);hh(1:end-1)]);  % put grey at bottom

function add_line(loc,clr)
   hold on
   t= loc(:,1); l = loc(:,2);
   plot(t,l,'-','Color',clr,'LineWidth',2)
   hold off
