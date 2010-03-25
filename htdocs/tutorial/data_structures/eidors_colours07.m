% Show lung images $Id$
img = inv_solve(imdl, zc_resp(:,1), zc_resp(:,20));
img.calc_colours.ref_level= 0;
img.calc_colours.cb_shrink_move = [0.5,0.8,.02];

clf;
subplot(221);

img.calc_colours.greylev = 0.01;
show_fem(img,[1,1]);
axis equal; axis off; print -r75 -dpng eidors_colours07a.png

img.calc_colours.greylev =  0.3;
show_fem(img,[1,1]);
axis equal; axis off; print -r75 -dpng eidors_colours07b.png

img.calc_colours.greylev = -0.01;
show_fem(img,[1,1]);
axis equal; axis off; print -r75 -dpng eidors_colours07c.png

img.calc_colours.cmap_type = 'draeger';
show_fem(img,[1,1]);
axis equal; axis off; print -r75 -dpng eidors_colours07d.png
