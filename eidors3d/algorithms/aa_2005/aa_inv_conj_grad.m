function img= aa_inv_conj_grad( inv_model, data1, data2)
% AA_INV_CONJ_GRAD inverse solver based on the CG
% inverse [Ref Shewchuck, 1994]
% img= aa_inv_conj_grad( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%

% $Id: aa_inv_conj_grad.m,v 1.1 2005-09-15 04:55:17 aadler Exp $

fwd_model= inv_model.fwd_model;
pp= aa_fwd_parameters( fwd_model );

% Note: no caching needed

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
sol = cg_inv( J'*W*J +  hp*R, J'*W*dva, imax, etol );

% create a data structure to return
img.name= 'solved by aa_inv_conj_grad';
img.elem_data = sol;
img.inv_model= inv_model;
img.fwd_model= fwd_model;

% CG code from [Shewchuck, 1994] Appendix B2
% www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf

function x= cg_inv( A, b, imax, etol )
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
   end
