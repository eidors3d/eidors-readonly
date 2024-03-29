function img= inv_solve_trunc_iterative( inv_model, data1, data2)
% INV_SOLVE_TRUNC_ITERATIVE using Morozov truncated iteration
% img= inv_solve_trunc_iterative( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% (C) 2005 David Stephenson. License: GPL version 2 or version 3
% $Id$

fwd_model= inv_model.fwd_model;

img_bkgnd= calc_jacobian_bkgnd( inv_model );
J = calc_jacobian(img_bkgnd);

% The one_step reconstruction matrix is cached
JtJ = eidors_cache(@calc_hessian, J, copt);

l_data1= length(data1); l1_0 = l_data1 ~=0;
l_data2= length(data2); l2_0 = l_data2 ~=0;
l_data= max( l_data1, l_data2 );


dva = calc_difference_data( data1, data2, fwd_model);

tol= 1e-4;
maxiter= 50;
if isfield(inv_model,'parameters');
    tol=     inv_model.parameters.term_tolerance;
    maxiter= inv_model.parameters.max_iterations;
end
sol = pcg(JtJ, J'*dva, tol, maxiter);


% create a data structure to return
img.name= 'solved by inv_solve_trunc_iterative';
img.elem_data = sol;
img.fwd_model= fwd_model;

function JtJ = calc_hessian(J)
   JtJ = J'*J;
