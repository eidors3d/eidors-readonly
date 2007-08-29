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
 
% REFS: 
% Andrea Borsic, Ph.D. Thesis, Oxford Brookes, UK 
% See Section 7.2.4 Lagged Diffusivity
% http://www.sc-aip.com/thesis/Andrea%20Borsic%20PhD.pdf
%
% W. Clem Karl, "Regularization in Image Restoration
% and Reconstruction," in Handbook of Image and Video
% Processing, A. Bovic, Ed. San Diego, CA: Academic Press,
% 2000, ch. 3.6, pp. 141-160.

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: aa_inv_total_var.m,v 1.10 2007-08-29 09:23:47 aadler Exp $

fwd_model= inv_model.fwd_model;
pp= aa_fwd_parameters( fwd_model );

img_bkgnd= calc_jacobian_bkgnd( inv_model );
J = calc_jacobian( fwd_model, img_bkgnd);

R = calc_R_prior( inv_model );
W = calc_meas_icov( inv_model );
hp= calc_hyperparameter( inv_model );


l_data1= length(data1); l1_0 = l_data1 ~=0;
l_data2= length(data2); l2_0 = l_data2 ~=0;
l_data= max( l_data1, l_data2 );

dva = calc_difference_data( data1, data2, inv_model.fwd_model);

imax= 3 ; % we need to do this avoid matlab version issues
if isfield(inv_model, 'parameters')
   if isfield(inv_model.parameters, 'max_iterations')
      imax= inv_model.parameters.max_iterations;
   end
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
   x= (H'*W*H + R'*R)\H'*W*y; %take reasonable first step
   Beta = 1e-4;
   for iter= 2:imax
      eidors_msg('AA_INV_TV: iteration %d',iter,3);
      WB = spdiags( 0.5 ./ sqrt( (R*x).^2 + Beta ), 0, n, n );
      x = (H'*W*H + R'*WB*R )\H'*W*y;
   end
