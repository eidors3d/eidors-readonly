% 2D solver $Id$

% Create a new inverse model, and set
% reconstruction model and fwd_model
imdl= mk_common_model('c2c2',16);
imdl.rec_model= cmdl;
fmdl.jacobian = @aa_calc_jacobian;
imdl.fwd_model= fmdl;

c2f= mk_coarse_fine_mapping( fmdl, cmdl);
imdl.fwd_model.coarse2fine = c2f;
imdl.RtR_prior = @gaussian_HPF_prior;
imdl.solve = @aa_inv_solve;
imdl.hyperparameter.value= 0.0001;


imgs= inv_solve(imdl, vh, vi);

show_fem(fmdl); ax= axis;
hold on
show_fem(imgs);
hold off
axis(ax); axis image
print_convert dm_geophys04a.png '-density 125'
