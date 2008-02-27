% Perturbation Jacobians
% $Id: perturb_jacobian_test.m,v 1.3 2008-02-27 02:06:39 aadler Exp $

% imdl= mk_common_model('c2c2',16);
  imdl= mk_common_model('a3cr',16);
  imdl= mk_common_model('n3r2',16);
  imdl.fwd_model.nodes = imdl.fwd_model.nodes*.25;
  img= calc_jacobian_bkgnd(imdl);

  img.fwd_model.normalize_measurements= 0;

  img.fwd_model.jacobian=   @np_calc_jacobian;
  img.fwd_model.system_mat= @np_calc_system_mat;
  img.fwd_model.solve=      @np_fwd_solve;
  J_np= calc_jacobian( img );

  img.fwd_model.jacobian=   @perturb_jacobian;
  img.fwd_model.system_mat= @np_calc_system_mat;
  img.fwd_model.solve=      @np_fwd_solve;
  J_np_p= calc_jacobian( img );

  img.fwd_model.jacobian=   @aa_calc_jacobian;
  img.fwd_model.system_mat= @aa_calc_system_mat;
  img.fwd_model.solve=      @aa_fwd_solve;
  J_aa= 2*calc_jacobian( img ); % 2 for bug in my code

  img.fwd_model.jacobian=   @perturb_jacobian;
  img.fwd_model.system_mat= @aa_calc_system_mat;
  img.fwd_model.solve=      @aa_fwd_solve;
  J_aa_p= 2*calc_jacobian( img ); % 2 for bug in my code

  norm(J_aa - J_aa_p,'fro')/norm(J_aa,'fro')
  norm(J_np - J_np_p,'fro')/norm(J_np,'fro')
  norm(J_np - J_aa_p,'fro')/norm(J_np,'fro')

