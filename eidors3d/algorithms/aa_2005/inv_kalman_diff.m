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
% $Id: inv_kalman_diff.m,v 1.4 2005-12-05 22:12:11 aadler Exp $

fwd_model= inv_model.fwd_model;
pp= aa_fwd_parameters( fwd_model );

% calc jacobian with homogeneous background
homg_img= eidors_obj('image', 'homog image', ...
                     'elem_data', ones( pp.n_elem ,1), ...
                     'fwd_model', fwd_model );

J = calc_jacobian( fwd_model, homg_img);

RtR = calc_RtR_prior( inv_model );
Q   = calc_meas_icov( inv_model );
hp  = calc_hyperparameter( inv_model );


l_data1= length(data1); l1_0 = l_data1 ~=0;
l_data2= length(data2); l2_0 = l_data2 ~=0;
l_data= max( l_data1, l_data2 );

dva= zeros(pp.n_meas, l_data);

if pp.normalize
   dva= 1 - data2 ./ data1;
else   
   dva= data1 - data2;
end

sol = kalman_inv( J, Q, hp^2*RtR, dva );

% create a data structure to return
img.name= 'solved by inv_kalman_diff';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;

% Kalman filter
% J - Jacobian NxM
% RegM - Regularization on the Measurements
% RegI - Regularization on the Image
function x= kalman_inv( J, RegM, RegI, y);
%x = (J'*RegM*J + RegI )\J'*RegM*y; return;
%Notation x_k1_k is x_{k+1|k}
 

% n is nmeas, m is ndata
[m,n]=  size(J);
% H is augmented matrix [J(x_k|k-1); RegI]
H_k1= [J;RegI];
% G (Gamma) is blockdiag [RegM, I]
scaling_const=1; % FIXME: what is this const?
G= speye(n+m)*scaling_const;
G(1:m,1:m) = RegM;
% F is the state transition matrix (I for random walk)
F_k= speye(n);
% Q is state noise covariance (model with I)
Q_k= speye(n);

% Initial C estimate. It is unitless, so no scale factor
C_k1_k1= speye(n);

% mean x_priori image - assume 0
x0= zeros(n,1);
RegI_x0= RegI*x0;
x_k1_k1= x0;

ll= size(y,2);
x= zeros(n,ll);
for i=1:ll
   eidors_msg('iteration %d',i,2);
   % Update variables
   C_k_k= C_k1_k1;
   x_k_k= x_k1_k1;
   % yi is [y_k; RegI*x0];
   yi= [y(:,i); RegI_x0];

   % Prediction
   x_k1_k = F_k * x_k_k;
   C_k1_k = F_k * C_k_k * F_k' + Q_k;
   % Correction
   HCHt   = H_k1 * C_k1_k * H_k1';
   K_k1   = C_k1_k * H_k1' / (HCHt + G);
   yerr   = yi - H_k1 * x_k1_k;
   x_k1_k1= x_k1_k + K_k1 * yerr; 
   C_k1_k1= (speye(n) - K_k1 * H_k1) * C_k1_k;

   % Store output
   x(:,i) = x_k1_k1;
end
   

