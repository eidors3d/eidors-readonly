function img= aa_inv_conj_grad( inv_model, data1, data2)
% AA_INV_CONJ_GRAD inverse solver based on the CG
% inverse [Ref Shewchuck, 1994]
% img= aa_inv_conj_grad( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: aa_inv_conj_grad.m,v 1.9 2007-08-29 09:20:53 aadler Exp $

fwd_model= inv_model.fwd_model;
pp= aa_fwd_parameters( fwd_model );

img_bkgnd= calc_jacobian_bkgnd( inv_model );
J = calc_jacobian( fwd_model, img_bkgnd);

R = calc_RtR_prior( inv_model );
W = calc_meas_icov( inv_model );
hp= calc_hyperparameter( inv_model );

maxiter= 50;
tol= 1e-4;
if isfield(inv_model,'parameters')
    tol =     inv_model.parameters.term_tolerance;
    maxiter = inv_model.parameters.max_iterations;
end



l_data1= length(data1); l1_0 = l_data1 ~=0;
l_data2= length(data2); l2_0 = l_data2 ~=0;
l_data= max( l_data1, l_data2 );

dva= zeros(pp.n_meas, l_data);

if pp.normalize
   dva= data2 ./ data1 - 1;
else   
   dva= data2 - data1;
end

n_img= size(dva,2);
sol = zeros( size(J,2), n_img );
Rx0 = zeros( size(R,1), 1);
for i=1:n_img
   sol(:,i) = cg_ls_inv0( J,  hp*R, dva(:,i), Rx0, maxiter, tol );
end

% create a data structure to return
img.name= 'solved by aa_inv_conj_grad';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;

% x = [J;R]\[y;R*x0] using Moore - Penrose inverse
% For comparison purposes
function x= cg_ls_inv1( J, R, y, Rx0, maxiter, tol )
   x = [J;R]\[y;Rx0];

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

   
