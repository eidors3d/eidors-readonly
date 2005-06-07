function img= np_inv_solve( inv_model, data1, data2)
% NP_INV_SOLVE inverse solver for Nick Polydorides EIDORS3D code
% img= np_inv_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time

% $Id: np_inv_solve.m,v 1.13 2005-06-07 03:01:42 aadler Exp $

% The one_step reconstruction matrix is cached
one_step_inv = eidors_obj('get-cache', inv_model, 'np_2003_one_step_inv');
if ~isempty(one_step_inv)
    eidors_msg('np_inv_solve: using cached value', 2);
else
    fwd_model= inv_model.fwd_model;
    p= np_fwd_parameters( fwd_model );

    % calc jacobian with homogeneous background
    homg_img= eidors_obj('image', 'homog image from np_inv_solve', ...
                         'elem_data', ones( p.n_elem ,1), ...
                         'fwd_model', fwd_model );

    J = calc_jacobian( fwd_model, homg_img);

    Reg = calc_image_prior( inv_model );
    tfac= calc_hyperparameter( inv_model );

    % Calculating a linear inverse solution
    one_step_inv= (J'*J +  tfac*Reg'*Reg)\J';

    eidors_obj('set-cache', inv_model, 'np_2003_one_step_inv', one_step_inv);
    eidors_msg('np_inv_solve: setting cached value', 2);
end


dva= data1.meas - data2.meas;
sol = one_step_inv * dva;

% create a data structure to return
img.name= 'solved by np_inv_solve';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;
