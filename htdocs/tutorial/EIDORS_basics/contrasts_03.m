   vh = fwd_solve(img);
   imgc = img;
   imgc.fwd_model.mdl_slice_mapper.npx = 128;
   imgc.fwd_model.mdl_slice_mapper.npy = 200;
   imgc.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
   imgc.calc_colours.ref_level = 1;
   q = show_current(imgc,vh.volt);
   hh=show_fem(imgc);
   set(hh,'EdgeColor',[1,1,1]*.75);
   hold on;

   sy = linspace(-2,2,20); sx= 0*sy - sz;
   hh=streamline(q.xp,q.yp, q.xc, q.yc,-sx,sy); set(hh,'Linewidth',2);

   hold off;

if ~exist('img_idx'); img_idx = 'a'; end
print_convert(sprintf('contrasts_03%c.png',img_idx));
img_idx = img_idx + 1;
