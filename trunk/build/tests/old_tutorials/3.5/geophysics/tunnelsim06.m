% Create 3D model of a tunnel $Id$ 

imdl = select_imdl( fmdl, {'Basic GN dif'});
imdl.rec_model = cmdl;
imdl.jacobian_bkgnd.value = cond_mdl;

% Do coarse2fine mapping. Rotate mdl to z dirn
f1mdl = fmdl; f1mdl.nodes = f1mdl.nodes(:,[2,3,1]);
f1mdl.mk_coarse_fine_mapping.z_depth = 1;
c2f= mk_coarse_fine_mapping( f1mdl, cmdl);
imdl.fwd_model.coarse2fine = c2f;

imdl.hyperparameter.value = 0.1;

imgr = inv_solve( imdl, vs_h, vs_i );

imgr.calc_colours.npoints= 128; subplot(221);
show_slices(imgr); print_convert tunnelsim06a.png
show_fem(imgr);    print_convert tunnelsim06b.png
show_fem(imgr);    axis(3*[-1,1,-1,1]); print_convert tunnelsim06c.png
