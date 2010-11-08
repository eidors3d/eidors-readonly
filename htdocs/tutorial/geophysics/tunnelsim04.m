% Create 3D model of a tunnel $Id$ 

% Reconstruct to a slice via coarse2fine
% Use a simple circular model without hole.
% Note that this will give the wrong electrode posns
imdl = mk_common_model('d2c2',N_elec);
imdl.rec_model = imdl.fwd_model;
imdl.rec_model.nodes = imdl.rec_model.nodes*5; % Enlarge
imdl.fwd_model = fmdl;
imdl.jacobian_bkgnd.value = cond_mdl;

% Do coarse2fine mapping. Rotate mdl to z dirn
f1mdl = fmdl; f1mdl.nodes = f1mdl.nodes(:,[2,3,1]);
f1mdl.mk_coarse_fine_mapping.z_depth = 1;
c2f= mk_coarse_fine_mapping( f1mdl, imdl.rec_model);
imdl.fwd_model.coarse2fine = c2f;

imdl.hyperparameter.value = 0.1;

imgr = inv_solve( imdl, vs_h, vs_i );

imgr.calc_colours.npoints= 128; subplot(221);
show_slices(imgr); print_convert tunnelsim04a.png
show_fem(imgr);    print_convert tunnelsim04b.png
