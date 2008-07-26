% Perturbation Jacobians
% $Id$

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
  J_aa= calc_jacobian( img ); % 2 for bug in my code

  img.fwd_model.jacobian=   @perturb_jacobian;
  img.fwd_model.system_mat= @aa_calc_system_mat;
  img.fwd_model.solve=      @aa_fwd_solve;
  J_aa_p= calc_jacobian( img ); % 2 for bug in my code

  norm(J_aa - J_aa_p,'fro')/norm(J_aa,'fro')
  norm(J_np - J_np_p,'fro')/norm(J_np,'fro')
  norm(J_np - J_aa_p,'fro')/norm(J_np,'fro')
  norm(J_np - J_aa  ,'fro')/norm(J_np,'fro')

% Demo model with EIDORS3D
  imdl= mk_common_model('n3r2',16);
  img= calc_jacobian_bkgnd(imdl);
  for i=1:16; imdl.fwd_model.electrode(i).z_contact= .01;end

% Calculate the normal Jacobian
  img.fwd_model.jacobian=   @np_calc_jacobian;
  img.fwd_model.system_mat= @np_calc_system_mat;
  img.fwd_model.solve=      @np_fwd_solve;
  J_np= calc_jacobian( img );

% Calculate the perturbation Jacobian
  vh= fwd_solve(img);

  for el= 1:50:size(img.fwd_model.elems,1);
     fprintf('\nelem#%03d: ',el);
     for delta= [1e-4,1e-6,1e-8] % delta is perturb
        img_i = img;
        img_i.elem_data(el)= img_i.elem_data(el) + delta;
        vi= fwd_solve(img_i);
        J_p = (vi.meas - vh.meas)/delta; % perturb Jacobian
        fprintf(' %8.6f', norm(J_p - J_np(:,el))/norm(J_np(:,el)));
     end
  end

