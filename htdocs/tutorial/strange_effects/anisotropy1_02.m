
% deformation
th = fmdl.nodes(:,1)/yt*(pi);
y = fmdl.nodes(:,2); y = y.*(th<=0) - y.*(th>0);
th = (th-pi/2).*(th<=0) + (pi/2-th).*(th>0);
[x,y] = pol2cart(th, y+4);
y = y+8.*(fmdl.nodes(:,1)<=0);
img2 = img; img2.fwd_model.nodes = [2*x, y];

show_fem(img2); axis image
img2.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img2);

img2.fwd_model.mdl_slice_mapper.npx = 200;
img2.fwd_model.mdl_slice_mapper.npy = 100;
q=  show_current(img2,vh.volt);
hh=show_fem(img2);
set(hh,'EdgeColor',[1,1,1]*.75);
hold on;
sx =  linspace(1.1,6.9,15); sy =  0*sx;
sy = -linspace(1.1,6.9,15); sx =  0*sy;
hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
hold off;
