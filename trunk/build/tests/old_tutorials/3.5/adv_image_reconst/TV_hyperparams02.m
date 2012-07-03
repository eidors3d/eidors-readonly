% TV: Reconstruction model $Id$

maxit=20;  % max number of iterations
imdl=mk_common_model('b2c0',16);

invtv= eidors_obj('inv_model', 'EIT inverse');
invtv.reconst_type= 'difference';
invtv.jacobian_bkgnd.value= 1;

invtv.fwd_model=                  imdl.fwd_model;
invtv.solve=                      @ab_tv_diff_solve;
invtv.R_prior=                    @ab_calc_tv_prior;
invtv.parameters.term_tolerance=  1e-6;
invtv.parameters.keep_iterations= 1;
invtv.parameters.max_iterations=  maxit;

subplot(221)
show_fem(invtv.fwd_model);
axis equal; axis off;
print_convert TV_hyperparams02a.png '-density 100';


