function img= aa_inv_solve( inv_model, data1, data2)
% AA_INV_SOLVE inverse solver for Nick Polydorides EIDORS3D code
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

% $Id: aa_inv_solve.m,v 1.11 2005-10-18 15:25:23 aadler Exp $

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

    J = calc_jacobian( fwd_model, homg_img );

    R = calc_image_prior( inv_model );
    W = calc_data_prior( inv_model );
    hp= calc_hyperparameter( inv_model );

    one_step_inv= (J'*W*J +  hp*R)\J'*W;

    eidors_obj('set-cache', inv_model, 'one_step_inv', one_step_inv);
    eidors_msg('aa_inv_solve: setting cached value', 2);
end

l_data1= length(data1); l1_0 = l_data1 ~=0;
l_data2= length(data2); l2_0 = l_data2 ~=0;
l_data= max( l_data1, l_data2 );

dva= zeros(pp.n_meas, l_data);

if pp.normalize
   dva= 1 - data2 ./ data1;
else   
   dva= data1 - data2;
end

sol = one_step_inv * dva;

% create a data structure to return
img.name= 'solved by aa_inv_solve';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;
