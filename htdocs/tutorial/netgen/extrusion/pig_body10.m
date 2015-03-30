cmdl.mk_coarse_fine_mapping.f2c_offset = [0,0,1];
cmdl.mk_coarse_fine_mapping.f2c_project = speye(3); % Scaling not required
cmdl.mk_coarse_fine_mapping.z_depth = 0.2;
c2f= mk_coarse_fine_mapping( fmdl, fmdlr);


imdl.name = 'CT pig 3D model';
imdl.fwd_model = fmdl;
imdl.rec_model = fmdlr;
imdl.fwd_model.coarse2fine = c2f;
imdl.jacobian_bkgnd.value = ones(size(fmdl.elems,1),1);
imdl.jacobian_bkgnd.value( fmdl.mat_idx{2} ) = 0.3;

imdl.fwd_model = mdl_normalize(imdl.fwd_model, 1);
imdl.hyperparameter.value = .03;
% Model background conductivity as lung
imr= inv_solve(imdl, vh, vi);

show_fem(imr); axis off ; axis tight

print_convert pig_body10a.jpg
