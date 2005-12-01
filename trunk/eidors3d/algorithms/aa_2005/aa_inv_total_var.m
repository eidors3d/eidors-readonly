function img= aa_inv_conj_grad( inv_model, data1, data2)
% AA_INV_CONJ_GRAD inverse solver based on the CG
% inverse [Ref Shewchuck, 1994]
% img= aa_inv_conj_grad( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: aa_inv_total_var.m,v 1.1 2005-12-01 09:17:05 aadler Exp $

fwd_model= inv_model.fwd_model;
pp= aa_fwd_parameters( fwd_model );

% calc jacobian with homogeneous background
homg_img= eidors_obj('image', 'homog image', ...
                     'elem_data', ones( pp.n_elem ,1), ...
                     'fwd_model', fwd_model );

J = calc_jacobian( fwd_model, homg_img);

R = calc_image_prior( inv_model );
W = calc_data_prior( inv_model );
hp= calc_hyperparameter( inv_model );


l_data1= length(data1); l1_0 = l_data1 ~=0;
l_data2= length(data2); l2_0 = l_data2 ~=0;
l_data= max( l_data1, l_data2 );

dva= zeros(pp.n_meas, l_data);

if pp.normalize
   dva= 1 - data2 ./ data1;
else   
   dva= data1 - data2;
end

imax= 100; etol= 1e-3;
n_img= size(dva,2);
sol = zeros( size(J,2), n_img );
for i=1:n_img
    % FIXME - ignore W for now
   sol(:,i) = tv_inv( J, hp*R, dva(:,i), imax, etol );
end

% create a data structure to return
img.name= 'solved by aa_inv_conj_grad';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;

function x= tv_inv( J, R, dva, imax, etol )
   x= (J'*J + R)\J'*dva;
