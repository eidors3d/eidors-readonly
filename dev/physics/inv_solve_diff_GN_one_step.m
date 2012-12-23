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

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

% TODO:
% Test whether Wiener filter form or Tikhonov form are faster
%  Tikhonov: RM= (J'*W*J +  hp^2*RtR)\J'*W;
%  Wiener:   P= inv(RtR); V = inv(W); RM = P*J'/(J*P*J' + hp^2*V)

dv = calc_difference_data( data1, data2, inv_model.fwd_model);
sol = get_RM( inv_model ) * dv;

img = physics_data_mapper(calc_jacobian_bkgnd( inv_model ));
img.name= 'solved by inv_solve_diff_GN_one_step';
img.elem_data = sol;
img.fwd_model= inv_model.fwd_model;
img = physics_data_mapper(img,1);

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
