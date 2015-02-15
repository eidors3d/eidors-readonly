   imgc = img;
   imgc.fwd_model.mdl_slice_mapper.npx = 128;
   imgc.fwd_model.mdl_slice_mapper.npy = 200;
   imgc.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
   imgc.calc_colours.ref_level = 1;
   q = show_current(imgc,vv.volt);
   hh=show_fem(imgc);
   set(hh,'EdgeColor',[1,1,1]*.75);
   hold on;

   sy = linspace(-2,2,20); sx= 0*sy - sz;
   hh=streamline(q.xp,q.yp, q.xc, q.yc,-sx,sy); set(hh,'Linewidth',2);

   hold off;

if ~exist('img_name'); img_name = '03a'; end
print_convert(sprintf('contrasts_%s.png',img_name));
