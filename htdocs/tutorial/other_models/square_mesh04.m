% 2D solver $Id$

% Create a new inverse model, and set
% reconstruction model and fwd_model
imdl= mk_common_model('c2c2',16);
imdl.rec_model= cmdl;
imdl.fwd_model= fmdl;

c2f= mk_coarse_fine_mapping( fmdl, cmdl);
imdl.fwd_model.coarse2fine = c2f;
%imdl.RtR_prior = @gaussian_HPF_prior;
imdl.RtR_prior = @noser_image_prior;
imdl.prior_use_fwd_not_rec = 1;
imdl.noser_image_prior.exponent= 0.5;
imdl.solve = @aa_inv_solve;
imdl.hyperparameter.value= 0.003;

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
print -dpng -r125 square_mesh04a.png
