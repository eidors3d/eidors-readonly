load('imdl_2x16_s4.mat');

ref_files = dir('*ref*.eit');
vh = [];
for rf = {ref_files.name}
   vh = [ vh mean(real(eidors_readdata(rf{1},'LQ3')),2)];
end
vh = mean(vh,2);


for L = [0 35 75]
   for R = [0 50 75 100]
      fname = sprintf('pos_L%d_R%d.eit',L,R);
      vi = mean(real(eidors_readdata(fname, 'LQ3')),2);
      if L==0 && R==0
         rimg = inv_solve(imdl_2x16_s4,vh,vi);
      else
         tmp =  inv_solve(imdl_2x16_s4,vh,vi);
         rimg.elem_data = [rimg.elem_data tmp.elem_data];
      end
   end
end
   rimg.calc_colours.ref_level = 0;
%    rimg.calc_colours.clim = max(abs(rimg.elem_data(:)))/2;
   rimg.show_slices.img_cols = size(rimg.elem_data,2);
   %%
   clear lvl
   lvl(:,3) = fliplr(.15:.02:.31);
   lvl(:,1:2) = Inf;
   show_slices(rimg,lvl);
   
   popt.supersampling_factor = 1;
   popt.resolution = 300;
   print_convert('tank_recon.png',popt);
   %%
   clf
   img = mk_image(imdl_2x16_s4.fwd_model,1);   
   h = show_fem(img);
   set(h, 'linewidth',.1);
   hold on
   [X, Y, Z] = sphere(40);
   X = 0.023*X; Y = 0.023*Y; Z = 0.023*Z;
   for L = [0 35 75]
      for R = [0 50 75 100]
         [x,y] = pol2cart(2*pi/32,R/100 * (0.14-0.024));
         surf(x+X,y+Y,Z+0.226+L/1000,'facecolor','b','edgecolor','none','facelighting','gouraud')
      end
   end
   l = linspace(0,2*pi,101);
   X = 0.14*sin(l);
   Y = 0.14*cos(l);
   for H = .14:.02:.32
      plot3(X,Y,H*ones(size(X)),'k','linewidth',2);
   end

   popt.supersampling_factor = 2;
   popt.resolution = 150;
   print_convert('tank.png',popt);