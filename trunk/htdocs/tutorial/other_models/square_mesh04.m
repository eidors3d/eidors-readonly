% 2D solver $Id$

% Create a new inverse model, and set
% reconstruction model and fwd_model
imdl= mk_common_model('c2c2',16);
imdl.rec_model= cmdl;
imdl.fwd_model= fmdl;

c2f= mk_coarse_fine_mapping( fmdl, cmdl);
imdl.fwd_model.coarse2fine = c2f;
%imdl.RtR_prior = @prior_gaussian_HPF;
imdl.RtR_prior = @prior_noser;
imdl.prior_use_fwd_not_rec = 1;
imdl.prior_noser.exponent= 0.5;
imdl.solve = @inv_solve_diff_GN_one_step;
imdl.hyperparameter.value= 0.05;

imgs= inv_solve(imdl, vh, vi);

show_fem(fmdl); ax= axis;
hold on
show_fem(imgs);
hold on;
h= trimesh(cmdl.elems,cmdl.nodes(:,1),cmdl.nodes(:,2));
set(h,'Color',[0,0,1]);
hold off
hold off
axis(ax); axis image
print_convert square_mesh04a.png
