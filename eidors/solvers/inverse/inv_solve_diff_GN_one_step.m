function img= inv_solve_diff_GN_one_step( inv_model, data1, data2)
% INV_SOLVE_DIFF_GN_ONE_STEP inverse solver using approach of Adler&Guardo 1996
% img= inv_solve_diff_GN_one_step( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix
%
% By default, the correct scaling of the solution that best fits the data
% is not calculated, possibly resulting in high solution errors reported by
% inv_solve. 
% To calculate the correct scaling, specify
%     inv_model.inv_solve_diff_GN_one_step.calc_step_size = 1;
% To provide a pre-calculated scaling, specify
%     inv_model.inv_solve_diff_GN_one_step.calc_step_size = 0;
%     inv_model.inv_solve_diff_GN_one_step.step_size = 0.8;
% The search for correct step_size is performed using FMINBND. The default
% search interval is [1e-5 1e1]. You can modify it by specifying:
%     inv_model.inv_solve_diff_GN_one_step.bounds = [10 200];
% Additional options for FMINBD can be passed as:
%     inv_model.inv_solve_diff_GN_one_step.fminbnd.MaxIter = 10;
%
% The optimal step_size is returned in img.step_size.
%
% See also INV_SOLVE, CALC_SOLUTION_ERROR, FMINBND

% (C) 2005-2013 Andy Adler and Bartlomiej Grychtol. 
% License: GPL version 2 or version 3
% $Id$

% TODO:
% Test whether Wiener filter form or Tikhonov form are faster
%  Tikhonov: RM= (J'*W*J +  hp^2*RtR)\J'*W;
%  Wiener:   P= inv(RtR); V = inv(W); RM = P*J'/(J*P*J' + hp^2*V)

dv = calc_difference_data( data1, data2, inv_model.fwd_model);
sol = get_RM( inv_model ) * dv;



img = data_mapper(calc_jacobian_bkgnd( inv_model ));
img.name= 'solved by inv_solve_diff_GN_one_step';
img.elem_data = sol;
img.fwd_model= inv_model.fwd_model;
img = data_mapper(img,1);

img = scale_to_fit_data(img, inv_model, data1, data2);


function RM = get_RM( inv_model )
   % The one_step reconstruction matrix is cached
   RM = eidors_obj('get-cache', inv_model, 'inv_solve_diff_GN_one_step');
   if ~isempty(RM)
       eidors_msg('inv_solve_diff_GN_one_step: using cached value', 3);
       return;
   end

   img_bkgnd= calc_jacobian_bkgnd( inv_model );
   J = calc_jacobian( img_bkgnd);

   RtR = calc_RtR_prior( inv_model );
   W   = calc_meas_icov( inv_model );
   hp  = calc_hyperparameter( inv_model );

   RM= (J'*W*J +  hp^2*RtR)\J'*W;

   eidors_obj('set-cache', inv_model, 'inv_solve_diff_GN_one_step', RM);
   eidors_msg('inv_solve_diff_GN_one_step: setting cached value', 3);
   
   
function [img step_size] = scale_to_fit_data(img, inv_model, data1, data2)
   % find the step size to multiply sol by to best fit data
   step_size = 1;
   do_step   = false;
   % If calc_step_size, ignore specified step_size
   try do_step = inv_model.inv_solve_diff_GN_one_step.calc_step_size; end
   
   if do_step
      eidors_msg('inv_solve_diff_GN_one_step: Calculating optimal step size to fit data',2);
      % options for fminbnd
      try 
         opt = inv_model.inv_solve_diff_GN_one_step.fminbnd;
      catch
         opt.Display = 'iter';
      end
      % range for fminbnd
      try
         range = inv_model.inv_solve_diff_GN_one_step.bounds;
      catch
         range = [1e-5 1e1];
      end
      step_size = fminbnd(@(x) to_optimize(img,inv_model,data1,data2, x), ...
                           range(1), range(2), opt);
   else
      % if not calculating, check if step_size provided
      try
         step_size = inv_model.inv_solve_diff_GN_one_step.step_size;
      end
   end
   img.elem_data = img.elem_data * step_size;
   img.step_size = step_size;

function out = to_optimize(img, inv_model, data1, data2, x)
   img.elem_data = img.elem_data*x;
   out = calc_solution_error(img, inv_model, data1, data2);

