% 2D solver $Id: centre_slice02.m,v 1.2 2008-03-25 14:37:50 aadler Exp $

imdl = mk_common_model('n3r2',[16,2]);
imdl.fwd_model = ng_mdl_16x2_coarse;

% Create coarse model
imdl2d= mk_common_model('c2c2',16);
c_mdl= imdl2d.fwd_model;
% fine model is offset by 15 in z dir. and shrink by 1/15
c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,15];
c_mdl.mk_coarse_fine_mapping.f2c_project = (1/15)*speye(3);
c2f= mk_coarse_fine_mapping( ng_mdl_16x2_coarse, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imdl.RtR_prior = @noser_image_prior;
img= inv_solve(imdl, vh, vi);
