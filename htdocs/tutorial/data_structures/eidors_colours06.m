% Show lung images $Id$
load montreal_data_1995
imdl = mk_common_model('c2t2',16);
imdl.fwd_model.normalize_measurements = 1;
imdl.RtR_prior= @gaussian_HPF_prior;
imdl.hyperparameter.value = 0.2;
img = inv_solve(imdl, zc_resp(:,1), zc_resp(:,20));

img.calc_colours.ref_level= 0;
img.calc_colours.cb_shrink_move = [0.5,0.8,.02];

subplot(221);
show_fem(img,[1,1]);
axis equal; axis off;

subplot(222);
img.calc_colours.ref_level =  0.1;
show_fem(img,[1,1]);
axis equal; axis off;

print -r75 -dpng eidors_colours06.png
