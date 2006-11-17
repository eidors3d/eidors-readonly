function img= inv_kalman_diff( inv_model, data1, data2)
% INV_KALMAN_DIFF inverse solver for difference EIT
% img= inv_kalman_diff( inv_model, data1, data2)
%
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% inv_model.inv_kalman_diff.sequence(1..nframes)
%            => sequence(1..nframes).meas_no
%            => sequence(1..nframes).time_delta (currently unused)
%            => sequence indicates which data are provided at
%               each frame (in order, starting with frame #1).
%               if not provided, all data are available at each frame
 
% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: inv_kalman_diff.m,v 1.10 2006-11-17 00:40:03 aadler Exp $

fwd_model= inv_model.fwd_model;
pp= aa_fwd_parameters( fwd_model );

img_bkgnd= calc_jacobian_bkgnd( inv_model );
J = calc_jacobian( fwd_model, img_bkgnd);

RtR = calc_RtR_prior( inv_model );
Q   = calc_meas_icov( inv_model );
hp  = calc_hyperparameter( inv_model );

if isfield(inv_model,'inv_kalman_diff')
   sequence= inv_model.inv_kalman_diff.sequence;
else
   sequence.meas_no= 1:size(J,1); % all data all the time
end


if pp.normalize
   dva= data2 ./ data1 - 1;
else   
   dva= data2 - data1;
end

 sol = kalman_inv( J, Q, hp^2*RtR, dva, sequence );
%R= calc_R_prior(inv_model);
%sol = kalman_inv_cgls( J, Q, hp^2*R, dva );

% create a data structure to return
img.name= 'solved by inv_kalman_diff';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;

% Kalman filter
% J - Jacobian NxM
% RegM - Regularization on the Measurements
% RegI - Regularization on the Image
% y is vector of measurements
% seq is sequence vector
function x= kalman_inv( J, RegM, RegI, y, seq);
%x = (J'*RegM*J + RegI )\J'*RegM*y; return;
%Notation x_k1_k is x_{k+1|k}
 

% n is nmeas, m is ndata
[m,n]=  size(J);
% G (Gamma) is blockdiag [RegM, I]
scaling_const=1; % FIXME: what is this const?

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
   eidors_msg('inv_kalman_diff: iteration %d',i,2);

   % H is augmented matrix [J(x_k|k-1); RegI]
   % need to account for sequence
   seq_i= seq( rem(i-1,length(seq))+1 ).meas_no;
   Jframe= J(seq_i,:); 
   H_k1= [Jframe;RegI];
   mm= size(Jframe,1);

   G= speye(n+mm)*scaling_const;
   G(1:mm,1:mm) = RegM(seq_i,seq_i);

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
%   C_k1_k1=  C_k1_k; - what is the effect of no update?

   % Store output
   x(:,i) = x_k1_k1;
end
   

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

   
