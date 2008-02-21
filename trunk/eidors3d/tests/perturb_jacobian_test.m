% Perturbation Jacobians
% $Id: perturb_jacobian_test.m,v 1.2 2008-02-21 21:15:21 aadler Exp $

  imdl= mk_common_model('a3cr',16);
  imdl.fwd_model.nodes = imdl.fwd_model.nodes*.25;
  img= eidors_obj('image','','elem_data',ones(768,1), ...
                  'fwd_model', imdl.fwd_model);

  img.fwd_model.normalize_measurements= 0;

  img.fwd_model.jacobian=   @np_calc_jacobian;
  img.fwd_model.system_mat= @np_calc_system_mat;
  img.fwd_model.solve=      @np_fwd_solve;
  J_np= calc_jacobian( img );

  img.fwd_model.jacobian=   @aa_calc_jacobian;
  img.fwd_model.system_mat= @aa_calc_system_mat;
  img.fwd_model.solve=      @aa_fwd_solve;
  J_aa= 2*calc_jacobian( img ); % 2 for bug in my code

  img.fwd_model.jacobian=   @perturb_jacobian;
  img.fwd_model.system_mat= @aa_calc_system_mat;
  img.fwd_model.solve=      @aa_fwd_solve;
  J_aa_p= 2*calc_jacobian( img ); % 2 for bug in my code

  norm(J_aa - J_aa_p,'fro')/norm(J_aa,'fro')
  norm(J_np - J_aa_p,'fro')/norm(J_np,'fro')

