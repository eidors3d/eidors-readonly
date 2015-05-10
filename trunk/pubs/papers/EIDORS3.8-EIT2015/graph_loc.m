function graph_loc
   %% Graph the lines of code
   loc
   clf
   add_line(loc_eidors, [2,0,0]/3);
   add_line(loc_htdocs, [0,2,0]/3);
   add_line(loc_other,  [0,0,2]/3);
   add_line(loc_dev,    [1,1,0]/3);
   k=1;for yy=2004:2015;
     tt{k} = (datenum(yy,1,1)-datenum(1970,1,1))*86400;
     tl{k} = num2str(yy);
   k=k+1; end
   set(gca,'xtick',cell2mat(tt))
   set(gca,'xticklabel',tl);
   ylim([0, ylim*[0;1]]);
   set(gca,'Box','off')
   set(gcf,'PaperPosition',[1,1,6,3]);
   %legend('Lung EIT','VILI/VALI','LPV','Location','NorthWest');
   %print -dpdf fig1_publications.pdf


function add_line(loc,clr)
   hold on
   t= loc(:,1); l = loc(:,2);
   plot(t,l,'-','Color',clr,'LineWidth',2)
   hold off
