% Perturbation Jacobians
% $Id: perturb_jacobian_test.m,v 1.1 2007-04-09 21:14:24 aadler Exp $

   imdl= mk_common_model('a3cr',16);
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
   J_aa= calc_jacobian( img );
   
   img.fwd_model.jacobian=   @perturb_jacobian;
   img.fwd_model.system_mat= @aa_calc_system_mat;
   img.fwd_model.solve=      @aa_fwd_solve;
   J_aa_p= calc_jacobian( img );
%  J_np_p= calc_jacobian( img );

