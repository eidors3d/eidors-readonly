% Calc Parameters $Id$

subplot(421); algs = get_list_of_algs;

for i= 1:length(algs)
   img = feval(algs{i}, vh, vi );
   param= [GREIT_noise_params(img, algs{i}, vh, vi); ... % noise parameters
           GREIT_sim_params(  img, xyzr_pt)];            % image parameters

   for j=1:size(param,1)
      plot(param(j,:)); set(gca,'XTickLabel',[])
      if     j==1; set(gca,'YLim',[0,1.5]);      % Noise Figure
      elseif j==2; set(gca,'YLim',[0,1.5]);      % Amplitude
      elseif j==3; set(gca,'YLim',[-0.05,0.25]);  % Posn Error
      elseif j==4; set(gca,'YLim',[0,0.4]);      % Resolution
      elseif j==5; set(gca,'YLim',[0,0.5]);      % Shape Deform
      elseif j==6; set(gca,'YLim',[0,1.5]);      % Ringing
      end
      print('-dpng','-r100',sprintf('simulation_test04_%d%d.png',i,j));
   end
end
