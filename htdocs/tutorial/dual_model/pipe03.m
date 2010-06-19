
% Create coarse model
imdl= mk_common_model('b2c2',16);
cmdl= imdl.fwd_model;

scl = 1;
cmdl.mk_coarse_fine_mapping.f2c_offset = [0,0,scl];
cmdl.mk_coarse_fine_mapping.f2c_project = (1/scl)*speye(3);
cmdl.mk_coarse_fine_mapping.z_depth = inf;
c2f= mk_coarse_fine_mapping( fmdl, cmdl);

% Create reconstruction model
imdl.rec_model= cmdl;
imdl.fwd_model= fmdl;
imdl.fwd_model.coarse2fine = c2f;
imdl.RtR_prior = @gaussian_HPF_prior;
imdl.solve = @aa_inv_solve;
imdl.hyperparameter.value= 1e-4;

imgc= inv_solve(imdl, vh, vi);

show_fem(imgc);
print_convert pipe03a.png '-density 75';
