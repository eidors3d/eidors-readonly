% Deformation $Id$
th = fmdl.nodes(:,1)/yt*(pi);
y = fmdl.nodes(:,2); y = y.*(th<=0) - y.*(th>0);
th = (th-pi/2).*(th<=0) + (pi/2-th).*(th>0);
[x,y] = pol2cart(th, y+4);
y = y+8.*(fmdl.nodes(:,1)<=0);
img2 = img; img2.fwd_model.nodes = [1.5*x, y];

img2.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img2);

img_v = rmfield(img2,'elem_data');
img_v.node_data = vh.volt;
hh=show_fem(img_v);
set(hh,'EdgeColor',[1,1,1]*.75);

q=  show_current(img2,vh.volt);
hold on;
sx =  linspace(1.1,6.9,15); sy =  0*sx;
sy = -linspace(1.1,6.9,15); sx =  0*sy;
hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
hold off;

axis image
print_convert anisotropy1_02a.png '-density 125'
