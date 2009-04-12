function img= aa_inv_solve( inv_model, data1, data2)
% AA_INV_SOLVE inverse solver using approach of Adler&Guardo 1996
% img= aa_inv_solve( inv_model, data1, data2)
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

fwd_model= inv_model.fwd_model;
pp= aa_fwd_parameters( fwd_model );

% The one_step reconstruction matrix is cached
one_step_inv = eidors_obj('get-cache', inv_model, 'aa_inv_solve');
if ~isempty(one_step_inv)
    eidors_msg('aa_inv_solve: using cached value', 3);
else
    img_bkgnd= calc_jacobian_bkgnd( inv_model );
    J = calc_jacobian( fwd_model, img_bkgnd);

    RtR = calc_RtR_prior( inv_model );
    W   = calc_meas_icov( inv_model );
    hp  = calc_hyperparameter( inv_model );

    one_step_inv= (J'*W*J +  hp^2*RtR)\J'*W;
%   We need to scale OSI to Vol*OSI*J*inv(diag(Vol)) = 1
%   Dvol = spdiags( pp.VOLUME, 0, pp.n_elem, pp.n_elem );
%   A = ( Dvol* one_step_inv) * (J/Dvol);
%   scl = A'\ones(pp.n_elem,1);
%   scl = (A*A'+ hp^2*RtR)\A*ones(pp.n_elem,1);
%   scl = (1:576)'; scl= 1 - (scl>200)*.4;
%   one_step_inv= spdiags(scl,0, pp.n_elem, pp.n_elem) * one_step_inv;

    eidors_obj('set-cache', inv_model, 'aa_inv_solve', one_step_inv);
    eidors_msg('aa_inv_solve: setting cached value', 3);
end

dv = calc_difference_data( data1, data2, inv_model.fwd_model);

sol = one_step_inv * dv;

% create a data structure to return
img.name= 'solved by aa_inv_solve';
img.elem_data = sol;
img.fwd_model= fwd_model;
