img.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img);

img.fwd_model.mdl_slice_mapper.npx = 200;
img.fwd_model.mdl_slice_mapper.npy = 100;
q = show_current(img,vh.volt);
hh=show_fem(img);
set(hh,'EdgeColor',[1,1,1]*.75);
hold on;
sy = linspace(-3,3,15); sx =  15+0*sy;
hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
hold off;
axis([-16,16,-3.3,3.3]);
save img2 img2; clear; eidors_cache clear; load img2 

img2.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img2);

img2.fwd_model.mdl_slice_mapper.npx = 200;
img2.fwd_model.mdl_slice_mapper.npy = 100;
q=  show_current(img2,vh.volt);
hh=show_fem(img2);
set(hh,'EdgeColor',[1,1,1]*.75);
hold on;
sx =  linspace(1.1,6.9,15); sy =  0*sx;
hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
hold off;
%axis([-16,16,-3.3,3.3]);
