function img= inv_kalman_diff( inv_model, data1, data2)
% INV_KALMAN_DIFF inverse solver for difference EIT
% img= inv_kalman_diff( inv_model, data1, data2)
%
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
 
% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: inv_kalman_diff.m,v 1.1 2005-12-05 13:14:50 aadler Exp $

fwd_model= inv_model.fwd_model;
pp= aa_fwd_parameters( fwd_model );

% calc jacobian with homogeneous background
homg_img= eidors_obj('image', 'homog image', ...
                     'elem_data', ones( pp.n_elem ,1), ...
                     'fwd_model', fwd_model );

J = calc_jacobian( fwd_model, homg_img);

RtR = calc_RtR_prior( inv_model );
Q = calc_meas_icov( inv_model );
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

n_img= size(dva,2);

sol = zeros( size(J,2), n_img );
for i=1:n_img
   sol(:,i) = kalman_inv( J, Q, hp*RtR, dva(:,i) );
end

% create a data structure to return
img.name= 'solved by inv_kalman_diff';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;

% FIXME: This is NOT a Kalman inv
function x= kalman_inv( H, Q, RtR, y);
   x = (H'*Q*H + RtR )\H'*Q*y;
