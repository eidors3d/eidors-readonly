function img= np_inv_solve( inv_model, data1, data2)
% NP_INV_SOLVE inverse solver for Nick Polydorides EIDORS3D code
% img= np_inv_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: np_inv_solve.m,v 1.29 2007-08-29 09:04:05 aadler Exp $

fwd_model= inv_model.fwd_model;

% The one_step reconstruction matrix is cached
one_step_inv = eidors_obj('get-cache', inv_model, 'np_2003_one_step_inv');
if ~isempty(one_step_inv)
    eidors_msg('np_inv_solve: using cached value', 2);
else
    p= np_fwd_parameters( fwd_model );

    img_bkgnd= calc_jacobian_bkgnd( inv_model );
    J = calc_jacobian( fwd_model, img_bkgnd);

    RtR = calc_RtR_prior( inv_model );
    hp= calc_hyperparameter( inv_model );

    % Calculating a linear inverse solution
    one_step_inv= (J'*J +  hp^2*RtR)\J';

    eidors_obj('set-cache', inv_model, 'np_2003_one_step_inv', one_step_inv);
    eidors_msg('np_inv_solve: setting cached value', 2);
end

if fwd_model.normalize_measurements
   dva= data2 ./ data1 - 1;
else   
   dva= data2 - data1;
end

sol = one_step_inv * dva;

% create a data structure to return
img.name= 'solved by np_inv_solve';
img.elem_data = sol;
img.fwd_model= fwd_model;
