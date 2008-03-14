% TV: Reconstruction model $Id: TV_hyperparams02.m,v 1.1 2008-03-14 15:57:32 aadler Exp $

maxit=40;  % max number of iterations
imdl=mk_common_model('b2c2',16);

invtv= eidors_obj('inv_model', 'EIT inverse');
invtv.reconst_type= 'difference';
invtv.jacobian_bkgnd.value= 1;

invtv.fwd_model=                  imdl.fwd_model;
invtv.solve=                      @ab_tv_diff_solve;
invtv.R_prior=                    @ab_calc_tv_prior;
invtv.parameters.term_tolerance=  1e-6;
invtv.parameters.keep_iterations= 1;
invtv.parameters.max_iterations=  maxit;


