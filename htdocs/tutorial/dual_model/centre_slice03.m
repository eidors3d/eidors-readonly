% 2D solver $Id: centre_slice03.m,v 1.3 2008-03-27 19:52:46 aadler Exp $


% Set coarse as reconstruction model
imdl.rec_model= c_mdl;
c_mdl.mk_coarse_fine_mapping.f2c_offset = [0,0,scl];
c_mdl.mk_coarse_fine_mapping.f2c_project = (1/scl)*speye(3);
c2f= mk_coarse_fine_mapping( f_mdl, c_mdl);
imdl.fwd_model.coarse2fine = c2f;
imdl.RtR_prior = @Gaussian_HPF_prior;
imdl.solve = @aa_inv_solve;
imdl.hyperparameter.value= 0.03;

imgc= inv_solve(imdl, vh, vi);

subplot(131)
show_slices(imgc);
