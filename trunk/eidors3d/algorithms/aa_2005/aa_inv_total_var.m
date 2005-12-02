function img= aa_inv_total_var( inv_model, data1, data2)
% AA_INV_TOTAL_VARIANCE inverse solver for difference EIT
% img= aa_inv_total_var( inv_model, data1, data2)
%
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% ALGORITHM: for image at iteration k x{k}
% Iterate:
%   x{k+1} = (H'*W*H + l^2*R'*WB( x{k} )*R )\H'*W*y
% For:
%   WB( x ) = diag( 0.5 ./ sqrt( (l*R*x).^2 + Beta ) )
% With Beta > 0 to approximate the L1 norm
%
% Comment on hyperparameter, we have l^2 in the iteration
% formula and l in WB, so overall hyperparameter is l. This
% makes the units work for the L1 norm.
 
% REF: 
% W. Clem Karl, "Regularization in Image Restoration
% and Reconstruction," in Handbook of Image and Video
% Processing, A. Bovic, Ed. San Diego, CA: Academic Press,
% 2000, ch. 3.6, pp. 141-160.

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: aa_inv_total_var.m,v 1.2 2005-12-02 09:11:21 aadler Exp $

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

if isfield(inv_model, 'parameters')
    imax= inv_model.parameters.max_iterations;
else
    imax= 3 ;
end
etol= 1e-3;
n_img= size(dva,2);

sol = zeros( size(J,2), n_img );
for i=1:n_img
   sol(:,i) = tv_inv( J, W, hp*R, dva(:,i), imax, etol );
end

% create a data structure to return
img.name= 'solved by aa_inv_conj_grad';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;

function x= tv_inv( H, W, R, y, imax, etol )
   n= size(R,1);
   x= zeros(n,1);
   Beta = 1e-8;
   for iter= 1:imax
      WB = spdiags( 0.5 ./ sqrt( (R*x).^2 + Beta ), 0, n, n );
      x = (H'*W*H + R'*WB*R )\H'*W*y;
   end
