% Image reconstruction of moving objects
% $Id: tutorial140a.m,v 1.1 2006-08-19 14:39:20 aadler Exp $

load netgen_moving_ball
vh= homg_tank;
vi= target_spirograph_fast;

imdl= mk_common_model('b2c',16);
imdl.hyperparameter.value= 1e-2;
imdl.RtR_prior= @laplace_image_prior;

% Use a Newton one-step solver
imdl.solve= @np_inv_solve;

subplot(211)
imgs1= inv_solve(imdl,vi,vh);
show_slices(imgs1(1:5),[ inf,inf,0,1,1])

% Use a Kalman filter solver
imdl.solve= @inv_kalman_diff;

subplot(212)
imgs2= inv_solve(imdl,vi,vh);
show_slices(imgs2(1:5),[ inf,inf,0,1,1])

% Print as images
pp=get(gcf,'paperposition');
set(gcf,'paperposition',pp.*[1,1,1,0.6]);
print -r100 -dpng tutorial140a.png;
set(gcf,'paperposition',pp.*[1,1,1,1]);
