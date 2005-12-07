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
% $Id: aa_inv_conj_grad.m,v 1.6 2005-12-07 22:45:04 aadler Exp $

fwd_model= inv_model.fwd_model;
pp= aa_fwd_parameters( fwd_model );

% Note: no caching needed

% calc jacobian with homogeneous background
homg_img= eidors_obj('image', 'homog image', ...
                     'elem_data', ones( pp.n_elem ,1), ...
                     'fwd_model', fwd_model );

J = calc_jacobian( fwd_model, homg_img);

R = calc_RtR_prior( inv_model );
W = calc_meas_icov( inv_model );
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
Rx0 = zeros( size(R,1), 1);
for i=1:n_img
%  sol(:,i) = cg_inv( J'*W*J +  hp^2*RtR, J'*W*dva(:,i), imax, etol );
   sol(:,i) = cg_ls_inv( chol(W)*J,  hp*R, dva(:,i), Rx0, imax, etol );
end

% create a data structure to return
img.name= 'solved by aa_inv_conj_grad';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;

% x = [J;R]\[y;R*x0] using Moore - Penrose inverse
function x= cg_ls_inv( J, R, y, Rx0, imax, etol )
  x = [J;R]\[y;Rx0]; % using Moore - Penrose inverse
  
   

% CG code from [Shewchuck, 1994] Appendix B2
% www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
function x= cg_inv_shewchuck( A, b, imax, etol )
   x= 1e-3*rand( size(A,2), 1);
   
   i=0;
   r= b- A*x;
   d= r;
   dnew= r'*r;
   d0= dnew;
   while (i<imax) & (dnew > etol^2*d0)
      q= A*d;
      a= dnew / (d'*q);
      x= x+ a*d;
      if rem(i,50)==0
          r= b- A*x;
      else
          r= r-a*q;
      end
      dold= dnew;
      dnew= r'*r;
      beta= dnew / dold;
      d= r+beta*d;
      i=i+1;
   end
