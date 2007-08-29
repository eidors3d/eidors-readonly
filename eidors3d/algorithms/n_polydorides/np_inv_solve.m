function img= np_inv_solve( inv_model, data1, data2)
% NP_INV_SOLVE inverse solver for Nick Polydorides EIDORS3D code
% img= np_inv_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
% inv_model.parameters.max_iterations (default 1);
% inv_model.parameters.term_tolerance (default 1e-3);

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: np_inv_solve.m,v 1.5 2007-08-29 09:15:32 aadler Exp $

[maxiter, tol] = get_parameters(inv_model);
  
dv = calc_difference_data( data1, data2, inv_model.fwd_model);

sol = one_step_inv_matrix(inv_model) * dv;

if maxiter>1
   RtR = calc_RtR_prior( inv_model );
   hp= calc_hyperparameter( inv_model );

   for iter=2:maxiter
      dv_sim= forward_solve_diff(inv_model, sol);
      eidors_msg('iter=%d, norm(err)= %f', iter, norm(dv_sim - dv),3);
      if norm(dv_sim - dv)<tol; break; end  % test tolerance

      img_bkgnd= calc_jacobian_bkgnd( inv_model );
      img_bkgnd.elem_data = img_bkgnd.elem_data + sol;
      J = calc_jacobian( inv_model.fwd_model, img_bkgnd);

      sol_upd= (J'*J +  hp^2*RtR)\(J' * (dv - dv_sim));
      sol = sol + sol_upd;
   end
end

% create a data structure to return
img.name= 'solved by np_inv_solve';
img.elem_data = sol;
img.fwd_model= inv_model.fwd_model;

function one_step_inv = one_step_inv_matrix(inv_model)
% The one_step reconstruction matrix is cached
   one_step_inv = eidors_obj('get-cache', inv_model, 'np_2003_one_step_inv');
   if ~isempty(one_step_inv)
       eidors_msg('np_inv_solve: using cached value', 2);
   else
       img_bkgnd= calc_jacobian_bkgnd( inv_model );
       J = calc_jacobian( inv_model.fwd_model, img_bkgnd);

       RtR = calc_RtR_prior( inv_model );
       hp= calc_hyperparameter( inv_model );

       % Calculating a linear inverse solution
       one_step_inv= (J'*J +  hp^2*RtR)\J';

       eidors_obj('set-cache', inv_model, 'np_2003_one_step_inv', one_step_inv);
       eidors_msg('np_inv_solve: setting cached value', 2);
   end

function dv_sim= forward_solve_diff(inv_model, sol) 
   img= calc_jacobian_bkgnd( inv_model );
   v_homg = fwd_solve(img);
   img.elem_data = img.elem_data + sol;
   v_solv = fwd_solve(img);
   dv_sim= calc_difference_data( v_homg, v_solv, inv_model.fwd_model);

function dva=diff_measurement(inv_model, data1, data2);
   fwd_model= inv_model.fwd_model;

   if fwd_model.normalize_measurements
      dva= data2 ./ data1 - 1;
   else   
      dva= data2 - data1;
   end

function [maxiter, tol] = get_parameters(inv_model);

   try
     maxiter= inv_model.parameters.max_iterations;
   catch
     maxiter= 1;
   end

   try
     tol = inv_model.parameters.term_tolerance;
   catch
     tol= 1e-3;
   end
