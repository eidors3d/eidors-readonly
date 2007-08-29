function img= inv_kalman_diff( inv_model, data1, data2)
% INV_KALMAN_DIFF inverse solver for difference EIT
% img= inv_kalman_diff( inv_model, data1, data2)
%
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% if inv_model.fwd_model.stimulation(:).delta_time
%   exists and is non_zero, then the kalman filter will
%   be applied to each data measurement separately
%
% Note that the classic Kalman filter assumes that the
%   time step between each measurement is constant
%   (ie as part of the state update eqn). inv_kalman_diff
%   cannot work with non-constant time steps
%
% if inv_model.inv_kalman_diff.keep_K_k1 = 1
%  then img outputs img.inv_kalman_diff.K_k1 = K_k1
%  this can be used to estimate noise properties
 
% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: inv_kalman_diff.m,v 1.1 2007-08-29 09:07:16 aadler Exp $

fwd_model= inv_model.fwd_model;
pp= aa_fwd_parameters( fwd_model );

img_bkgnd= calc_jacobian_bkgnd( inv_model );
J = calc_jacobian( fwd_model, img_bkgnd);

RtR = calc_RtR_prior( inv_model );
Q   = calc_meas_icov( inv_model );
hp  = calc_hyperparameter( inv_model );

if isfield(fwd_model.stimulation(1),'delta_time')
   delta_time= [fwd_model.stimulation(:).delta_time];
   if diff(delta_time) ~= 0;
      error('All time steps must be same for kalman filter');
   end
else
   delta_time=0;
end

% sequence is a vector location of each stimulation in the frame
if delta_time == 0
   sequence = size(J,1);
else
   for i=1:length(fwd_model.stimulation)
      sequence(i) = size(fwd_model.stimulation(i).meas_pattern,1);
   end
   sequence= cumsum( sequence );
end


if pp.normalize
   dva= data2 ./ data1 - 1;
else   
   dva= data2 - data1;
end

[sol, K_k1] = kalman_inv( J, hp^2*Q, RtR, dva, sequence );

% create a data structure to return
img.name= 'solved by inv_kalman_diff';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;

try % keep parameter if requested
   if inv_model.inv_kalman_diff.keep_K_k1
      img.inv_kalman_diff.K_k1 = K_k1;
   end
end

% Kalman filter - estimates x where z_k = H_k*x_k + noise
% J - Jacobian NxM
% RegM - Regularization on the Measurements
% RegI - Regularization on the Image
% y is vector of measurements
% seq is sequence vector
%
% K_k1 is the linearized reconstruction matrix for
%   the final step. It can be used to estimate noise
%   properties of the algorithm 
%
function [x, K_k1]= kalman_inv( J, RegM, RegI, y, seq);
%x = (J'*RegM*J + RegI )\J'*RegM*y; return;
%Notation x_k1_k is x_{k+1|k}

% n is nmeas, m is ndata
[m,n]=  size(J);

% F is the state transition matrix (I for random walk)
F_k= speye(n);
% Q is state noise covariance (model with I)
Q_k= RegI;

% Initial error covariance estimate.
C_k1_k1= speye(n);

% mean x_priori image - assume 0
x0= zeros(n,1);
x_k1_k1= x0;

ll= size(y,2);
x= zeros(n,ll*length(seq));

seq= [0;seq(:)];
iter=0;
for i=1:ll
   for ss= 2:length(seq);
      eidors_msg('inv_kalman_diff: iteration %d.%d',i,ss-1,2);

      seq_i= (seq(ss-1)+1) : seq(ss);

% The Kalman filter doesn't need the regularization at all
      H_k1= J(seq_i,:);
      yi= y(seq_i,i);
      G = RegM(seq_i,seq_i);
      [x_k1_k1, C_k1_k1, K_k1] = ...
             kalman_step( x_k1_k1, C_k1_k1, ...
                          H_k1, yi, F_k, Q_k, G );
      iter=iter+1;
      x(:,iter) = x_k1_k1;
   end
end

function [x_k1_k1, C_k1_k1, K_k1] = ...
                  kalman_step( x_k_k, C_k_k, ...
                               H_k1, yi, F_k, Q_k, G )
   n= size(H_k1,2);

   % Prediction
   x_k1_k = F_k * x_k_k;
   C_k1_k = F_k * C_k_k * F_k' + Q_k;
   % Correction
   HCHt   = H_k1 * C_k1_k * H_k1';
   K_k1   = C_k1_k * H_k1' / (HCHt + G);
   yerr   = yi - H_k1 * x_k1_k;
   x_k1_k1= x_k1_k + K_k1 * yerr; 
   C_k1_k1= (speye(n) - K_k1 * H_k1) * C_k1_k;
   

function x= kalman_inv_cgls( J, RegM, RegI, y);
   [m,n]=  size(J);
   Rx0= zeros(n,1);
   H= [chol(W)*J;RegI];
   b= [y;Rx0];
   x= H\b;

% Adapted from code in Hansen's regularization tools
function x= cg_ls_inv0( J, R, y, Rx0, maxiter, tol )
%  A = [J;R];
   b=[y;Rx0];
   [m,n]= size(J);
   m_idx= 1:m; n_idx = m+(1:n);
   Jt = J.';
   Rt = R.';
   x = zeros(n,1); 
%  d = A'*b; 
   d = Jt*b(m_idx) + Rt*b(n_idx);
   r = b; 
   normr2 = d'*d; 
    
   k=0; % Iterate. 
   x_delta_filt= 1; x_filt= .1;
   while 1 
     % Update x and r vectors. 
%    Ad = A*d;
     Ad = [J*d;R*d];
     Alpha = normr2/(Ad'*Ad); 
     xpre= x;
     x  = x + Alpha*d; 

     k= k+1; if k==maxiter; break ; end

     x_delta= norm(xpre-x)/norm(x);
     x_delta_filt= x_delta_filt*(1-x_filt) + x_filt*x_delta;

     if x_delta_filt<tol; break ; end

     r  = r - Alpha*Ad; 
%    s  = A'*r; 
     s  = Jt*r(m_idx) + Rt*r(n_idx);
    
     % Update d vector. 
     normr2_new = s'*s; 
     Beta = normr2_new/normr2; 
     normr2 = normr2_new; 
     d = s + Beta*d; 
      
   end 
%     disp([k, x_delta, x_delta_filt]);

   
