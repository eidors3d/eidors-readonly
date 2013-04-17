%Inverse solution
imdl = mk_common_model('b2c2',32); %generic mdl
imdl.solve = @inv_solve_abs_GN; %Default Gauss Newton solvers
imdl.fwd_model = mdl_i;
imdl.reconst_type = 'absolute';
imdl.jacobian_bkgnd.value= cond_h;

imdl.parameters.show_iterations=1; %Show iteration progress
imdl.parameters.max_iterations = 3 ; %Number of iterations

cb.cb_shrink_move = [0.5,0.8,0.00];

   imdl.RtR_prior=@prior_laplace; hp = 1e-4^2;
%  imdl.RtR_prior=@prior_noser;   hp = 1e-2^2;

imdl.hyperparameter.value = hp;
img   = inv_solve(imdl, v_i);
img_n = inv_solve(imdl, v_i_n);

clf; show_3d_slices(img,0.6,0.3,0.3);
eidors_colourbar(setfield(img,'calc_colours',cb));
print_convert compare_3D_abs_GN_02a.png '-density 75'

clf; show_3d_slices(img_n,0.6,0.3,0.3);
eidors_colourbar(setfield(img_n,'calc_colours',cb));
print_convert compare_3D_abs_GN_02b.png '-density 75'

%img.show_slices.levels = [inf,inf,.5,1,1;inf,inf,0.6,2,1]; show_slices(img);
