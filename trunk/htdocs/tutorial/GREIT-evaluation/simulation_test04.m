% Calc Parameters $Id$

subplot(421); algs = get_list_of_algs;

for i= 1:length(algs)
   img = feval(algs{i}, vh, vi );
   param= GREIT_sim_params( img, xyzr_pt);

   noise = 0.01*std(vh)*randn(208,1000);
   snr_y = mean(abs(vi-vh*ones(1,size(vi,2))),1) / mean(std(noise,[],1),2); 

   im_n= feval(algs{i}, vh, vh*ones(1,size(noise,2)) + noise);
   snr_x = mean(mean(abs(img),1),2) / mean(abs(im_n(:)));
   param= [param; [snr_x(:)./snr_y(:)]' ];

   for j=1:size(param,1)
      plot(param(j,:)); set(gca,'XTickLabel',[])
      if     j==1; set(gca,'YLim',[0,1.3]);
      elseif j==2; set(gca,'YLim',[-0.05,0.2]);
      elseif j==3; set(gca,'YLim',[0,0.4]);
      elseif j==4; set(gca,'YLim',[0,0.5]);
      elseif j==5; set(gca,'YLim',[0,1.0]);
      elseif j==6; set(gca,'YLim',[0,4.0]);
      end
      print('-dpng','-r100',sprintf('simulation_test04_%d%d.png',i,j));
   end
end
