% 2D solver $Id$

% Create a new inverse model, and set
% reconstruction model and fwd_model
imdl= mk_common_model('c2c2',16); im_fm= imdl.fwd_model;
imdl.rec_model= s_mdl;
imdl.fwd_model= f_mdl;
imdl.fwd_model.stimulation= im_fm.stimulation;
imdl.fwd_model.meas_select= im_fm.meas_select;

s_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,5];
s_mdl.mk_coarse_fine_mapping.f2c_project = (1/15)*speye(3);
s_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, s_mdl);
imdl.fwd_model.coarse2fine = c2f;
imdl.RtR_prior = @gaussian_HPF_prior;
imdl.solve = @aa_inv_solve;
imdl.hyperparameter.value= 0.001;

imgs= inv_solve(imdl, vh, vi);
