% Basic Image reconstruction
% $Id: tutorial110a.m,v 1.1 2006-08-19 04:09:40 aadler Exp $

% Load some data
load netgen_moving_ball

% Get a 2D image reconstruction model
imdl= mk_common_model('c2c');
vi= target_spiral(50); vh= homg_tank;

imdl.hyperparameter.value= 1e-3;
img= inv_solve(imdl, vi, vh);
subplot(231); show_slices(img);
subplot(234); mesh( show_slices(img) );
view(-11,44); axis([0 128 0 128 -.0002 .002]);

imdl.hyperparameter.value= 1e-2;
img= inv_solve(imdl, vi, vh);
subplot(232); show_slices(img);
subplot(235); mesh( show_slices(img) );
view(-11,44); axis([0 128 0 128 -.0002 .002]);

imdl.hyperparameter.value= 1e-1;
img= inv_solve(imdl, vi, vh);
subplot(233); show_slices(img);
subplot(236); mesh( show_slices(img) );
view(-11,44); axis([0 128 0 128 -.0002 .002]);


print -r100 -dpng tutorial110a.png;
