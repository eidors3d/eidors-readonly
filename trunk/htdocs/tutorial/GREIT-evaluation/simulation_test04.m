% Calc Parameters $Id$

subplot(421); algs = get_list_of_algs;

for i= 1:length(algs)
   img = feval(algs{i}, vh, vi );
   param= GREIT_sim_params( img, xyzr_pt);
   for j=1:size(param,1)
      plot(param(j,:)); set(gca,'XTickLabel',[])
      if     j==1; set(gca,'YLim',[0,1.2]);
      elseif j==2; set(gca,'YLim',[0,0.2]);
      elseif j==3; set(gca,'YLim',[0,0.4]);
      end
      print('-dpng','-r100',sprintf('simulation_test04_%d%d.png',i,j));
   end
end
