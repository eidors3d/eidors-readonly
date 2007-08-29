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
% $Id: inv_solve_trunc_iterative.m,v 1.2 2007-08-29 09:10:12 aadler Exp $

fwd_model= inv_model.fwd_model;

img_bkgnd= calc_jacobian_bkgnd( inv_model );
J = calc_jacobian( fwd_model, img_bkgnd);

% The one_step reconstruction matrix is cached
JtJ = eidors_obj('get-cache', inv_model, 'Hessian');
if ~isempty(JtJ)
    eidors_msg('inv_solve_trunc_iterative: using cached value', 2);
else
    JtJ= J'*J;

    eidors_obj('set-cache', inv_model, 'Hessian', JtJ);
    eidors_msg('inv_solve_trunc_iterative: setting cached value', 2);
end

l_data1= length(data1); l1_0 = l_data1 ~=0;
l_data2= length(data2); l2_0 = l_data2 ~=0;
l_data= max( l_data1, l_data2 );

dva= zeros(size(J,1), l_data);

if pp.normalize
   dva= data2 ./ data1 - 1;
else   
   dva= data2 - data1;
end

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
