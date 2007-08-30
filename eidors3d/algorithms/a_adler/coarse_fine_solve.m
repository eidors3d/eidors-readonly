function img= coarse_fine_solve( inv_model, data1, data2)
% COARSE_FINE_SOLVER
% img= np_inv_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% Inv_model_parameters
% inv_model.coarse_fine.solve   => actual underlying solver
% inv_model.coarse_fine.mapping => sparse matrix to map coarse->fine

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: coarse_fine_solve.m,v 1.9 2007-08-30 03:37:02 aadler Exp $

fwd_model= inv_model.fwd_model;

% The one_step reconstruction matrix is cached
one_step_inv = eidors_obj('get-cache', inv_model, 'coarse_fine_one_step_inv');
if ~isempty(one_step_inv)
    eidors_msg('coarse_fine_solve: using cached value', 2);
else

%TODO: Note that this is a work in progress
% The right way to do it is to modify the
% call to calc_jacobian and RtR_prior, so that
% they're modified by M, as shown here. Then the
% inverse can proceed normally.

    p= np_fwd_parameters( fwd_model );

    img_bkgnd= calc_jacobian_bkgnd( inv_model );
    J = calc_jacobian( fwd_model, img_bkgnd);

    RtR = calc_RtR_prior( inv_model );
    hp= calc_hyperparameter( inv_model );

    M= inv_model.coarse_fine.mapping;
    JM= J*M;
    MtRtRM= M'*RtR*M;

    % Calculating a linear inverse solution
    one_step_inv= (JM'*JM +  hp^2*MtRtRM)\JM';

    eidors_obj('set-cache', inv_model, 'coarse_fine_one_step_inv', one_step_inv);
    eidors_msg('coarse_fine_solve: setting cached value', 2);
end

dva= calc_difference_data(data1, data2, fwd_model);

sol = one_step_inv * dva;

% create a data structure to return
img.name= 'solved by coarse_fine_solve';
img.elem_data = inv_model.coarse_fine.mapping * sol;
img.fwd_model= fwd_model;

