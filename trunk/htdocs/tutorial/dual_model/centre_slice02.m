% 2D solver $Id: centre_slice02.m,v 1.3 2008-03-25 19:21:49 aadler Exp $

imdl = mk_common_model('b3cr',[16,2]); f_mdl=imdl.fwd_model;
%imdl.fwd_model = ng_mdl_16x2_coarse;

% Create coarse model
imdl2d= mk_common_model('b2c2',16);
c_mdl= imdl2d.fwd_model;
% if fine model is offset use:
%c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,15];
%c_mdl.mk_coarse_fine_mapping.f2c_project = (1/15)*speye(3);
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imdl.RtR_prior = @noser_image_prior;
imdl.solve = @aa_inv_solve;
imdl.hyperparameter.value= 0.1;

img= inv_solve(imdl, vh, vi);
imgc= img; imgc.fwd_model= c_mdl;

subplot(131)
show_slices(imgc);
