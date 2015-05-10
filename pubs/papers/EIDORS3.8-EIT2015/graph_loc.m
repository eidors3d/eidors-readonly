function graph_loc
   %% Graph the lines of code
   loc
   clf
if 0 % From python git-loc ... not sure we trust the data
   add_line(loc_eidors, [2,0,0]/3);
   add_line(loc_htdocs, [0,2,0]/3);
   add_line(loc_other,  [0,0,2]/3);
   add_line(loc_dev,    [1,1,0]/3);
end
   hold on;
   plot(dat(:,2),1e-3*dat(:,3),            '-','Color',[2,0,0]/3,'LineWidth',2);
   plot(dat(:,2),1e-3*dat(:,3:5)*[1;-1;-1],'-','Color',[1,1,0]/3,'LineWidth',2);
   plot(dat(:,2),1e-3*dat(:,4),            '-','Color',[0,2,0]/3,'LineWidth',2);
   plot(dat(:,2),1e-3*dat(:,5),            '-','Color',[0,0,2]/3,'LineWidth',2);

   k=1;for yy=2004:2016;
     tt(k) = (datenum(yy,1,1)-datenum(1970,1,1))*86400;
     tl{k} = num2str(yy); if rem(yy,2)==1; tl{k}=''; end
   k=k+1; end
   set(gca,'xtick',tt)
   set(gca,'xticklabel',tl);
   ylim([0, ylim*[0;1]]);
release= [
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
   plot([1;1]*release(:)',[0,ylim*[0;1]], 'Color',[2,2,2]/4);
   hold off;
   set(gca,'Box','off')
   set(gcf,'PaperPosition',[1,1,6,2.5]);
   legend('Total','Eidors','Tutorials','Dev','location','northwest');
   ylabel('Lines of Code (x1000)')
   print -dpdf fig_loc.pdf


function add_line(loc,clr)
   hold on
   t= loc(:,1); l = loc(:,2);
   plot(t,l,'-','Color',clr,'LineWidth',2)
   hold off
