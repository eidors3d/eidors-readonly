sz = 5; img_idx = 'b';
for ellipse_x = [0.5,1,2];
   img = contrasts_04_modeller( sz, ellipse_x); 
   targ = img.fwd_model.mat_idx{1};
   vr = fwd_solve(img);
   imgv= rmfield(img,'elem_data');
   for contrast = linspace( -2,2,5);
      img.elem_data( targ ) = exp(contrast);
      vv = fwd_solve(img);
      vv.volt = vv.volt - vr.volt;
%     contrasts_03;
      img.fwd_model.mdl_slice_mapper.npx = 128;
      img.fwd_model.mdl_slice_mapper.npy = 200;
      q = show_current(img,vv.volt);

   imgv.node_data = vv.volt;
   imgv.calc_colours.ref_level = mean(vv.volt);
   hh=show_fem(imgv);
   set(hh,'EdgeColor',[1,1,1]*.75);
   hold on;

   sy = linspace(-2,2,20); sy(abs(sy)<0.5/ellipse_x) = [];
   sx= 0*sy;
   hh=streamline(q.xp,q.yp, q.xc, q.yc, sx, sy); set(hh,'Linewidth',2);
   hh=streamline(q.xp,q.yp,-q.xc,-q.yc, sx, sy); set(hh,'Linewidth',2);

   hold off;
   img_name = sprintf('05%c',img_idx); img_idx= img_idx+1;
   print_convert(sprintf('contrasts_%s.png',img_name));
   end
end

