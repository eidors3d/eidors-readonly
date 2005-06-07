function img= aa_inv_solve( inv_model, data1, data2)
% NP_INV_SOLVE inverse solver for Nick Polydorides EIDORS3D code
% img= aa_inv_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time

% $Id: aa_inv_solve.m,v 1.2 2005-06-07 03:18:08 aadler Exp $

fwd_model= inv_model.fwd_model;
pp= aa_fwd_parameters( fwd_model );

% The one_step reconstruction matrix is cached
one_step_inv = eidors_obj('get-cache', inv_model, 'one_step_inv');
if ~isempty(one_step_inv)
    eidors_msg('aa_inv_solve: using cached value', 2);
else
    % calc jacobian with homogeneous background
    homg_img= eidors_obj('image', 'homog image', ...
                         'elem_data', ones( pp.n_elem ,1), ...
                         'fwd_model', fwd_model );

    J = calc_jacobian( fwd_model, homg_img);

    R = calc_image_prior( inv_model );
    hp= calc_hyperparameter( inv_model );

    one_step_inv= (J'*J +  hp*R)\J';

    eidors_obj('set-cache', inv_model, 'one_step_inv', one_step_inv);
    eidors_msg('aa_inv_solve: setting cached value', 2);
end

if pp.normalize
   dva= 1 - data2.meas ./ data1.meas;
else   
   dva= data1.meas - data2.meas;
end

sol = one_step_inv * dva(:);

% create a data structure to return
img.name= 'solved by aa_inv_solve';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;
