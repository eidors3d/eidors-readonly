% 2D solver $Id: centre_slice06.m,v 1.1 2008-03-27 20:56:01 aadler Exp $


% Set coarse as reconstruction model
imdl.rec_model= c_mdl;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imdl.RtR_prior = @gaussian_HPF_prior;
imdl.solve = @aa_inv_solve;
imdl.hyperparameter.value= 0.03;

imgc= inv_solve(imdl, vh, vi);
subplot(131); show_slices(imgc);

c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,-.3];
c_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imgc0= inv_solve(imdl, vh, vi);

imgc0= inv_solve(imdl, vh, vi);
subplot(132); show_slices(imgc0);

c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,.3];
c_mdl.mk_coarse_fine_mapping.z_depth = 0.1;
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;

imgc1= inv_solve(imdl, vh, vi);
subplot(133); show_slices(imgc1);

print -r150 -dpng centre_slice06a.png;
