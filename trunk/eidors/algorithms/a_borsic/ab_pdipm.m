function img=ab_pdipm( inv_model, data1, data2)
% AB_PDIPM inverse solver using Primal/Dual interior point method
% img= ab_pdipm( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
%  inv_model.ab_pdipm.norm_data  1 or 2 (DEFAULT 2)
%  inv_model.ab_pdipm.norm_prior 1 or 2 (DEFAULT 2)
%  inv_model.ab_pdipm.beta     (default 1e-6)
%
% Parameters:
%  max_iters =  inv_model.parameters.max_iteration (default 10)
%      Max number of iterations before stopping
%  min change = inv_model.parameters.min_change   (default 0)
%      Min Change in objective fcn (norm(y-Jx)^2 + hp*TV(x)) before stopping
% beta is the parameter that smooths the TV functional

% (C) 2008 Andrea Borsic. License: GPL version 2 or version 3
% $Id$


pp= process_parameters(inv_model);

fwd_model= inv_model.fwd_model;

d=calc_difference_data( data1, data2, fwd_model);

img_bkgnd=calc_jacobian_bkgnd( inv_model );
J=calc_jacobian( fwd_model, img_bkgnd);

alpha=calc_hyperparameter( inv_model );
L=calc_R_prior( inv_model );
W= calc_meas_icov( inv_model );


if     pp.norm_data==2 && pp.norm_image==2
  x= pdipm_2_2( J,W,alpha*L,d, pp);;
elseif pp.norm_data==2 && pp.norm_image==1
  x= pdipm_2_1( J,W,alpha*L,d, pp);;
elseif pp.norm_data==1 && pp.norm_image==2
  x= pdipm_1_2( J,W,alpha*L,d, pp);;
elseif pp.norm_data==1 && pp.norm_image==1
  x= pdipm_1_1( J,W,alpha*L,d, pp);;
end

% create a data structure to return
img.name = 'ab_pdipm';
img.elem_data = x;
img.fwd_model = fwd_model;

function s= pdipm_2_2( J,W,L,d, pp);
   s= (J'*W*J + L'*L)\J'*W*d;

function s= pdipm_1_2( J,W,L,d, pp);
   [M,N] = size(J); % M measurements, N parameters
   s= zeros( N, 1 ); % solution - start with zeros
   x= zeros( M, 1 ); % dual var - start with zeros

   for loop = 1:pp.max_iter
      % Define variables
      f = J*s - d;             F= spdiags(f,0,M,M);
                               X= spdiags(x,0,M,M);
      e = sqrt(f.^2 + pp.beta);E= spdiags(e,0,M,M);

      % Define derivatives
      dFc_ds = (speye(M,M) - X*inv(E)*F)*J;
      dFc_dx = -E;
      dFf_ds = L'*L;
      dFf_dx = J'*W;

      dsdx = -[dFc_ds, dFc_dx; dFf_ds, dFf_dx] \ ...
              [ f-E*x; J'*W*x + L'*L*s ];

      s= s + dsdx(1:N);
      x= x + dsdx(N+ (1:M));
      x( x> 1 ) = 1;
      x( x<-1 ) = -1;
fprintf('.\n');
   end

function pp= process_parameters(imdl);
   try    pp.max_iter = imdl.parameters.max_iterations;
   catch  pp.max_iter = 10;
   end

   try    pp.min_change = imdl.parameters.min_change;
   catch  pp.min_change = 0;
   end

   try    pp.beta = imdl.ab_pdipm.beta; 
   catch  pp.beta = 1e-6;
   end

   try    pp.norm_data = imdl.ab_pdipm.norm_data;
   catch  pp.norm_data = 2;
   end

   try    pp.norm_image = imdl.ab_pdipm.norm_image;
   catch  pp.norm_image = 2;
   end


