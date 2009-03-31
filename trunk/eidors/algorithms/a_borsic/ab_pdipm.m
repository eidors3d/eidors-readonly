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

      ds = dsdx(1:N);
      dx = x_update(x, dsdx(N+(1:M)));

      s= s + ds; x= x + dx;
fprintf('+');
   end

function s= pdipm_2_1( J,W,L,d, pp);
   [M,N] = size(J);   % M measurements, N parameters
   [G  ] = size(L,1); % E edges
   s= zeros( N, 1 ); % solution - start with zeros
   x= zeros( G, 1 ); % dual var - start with zeros

   for loop = 1:pp.max_iter
      % Define variables
      f = L*s;                 F= spdiags(f,0,G,G);
                               X= spdiags(x,0,G,G);
      e = sqrt(f.^2 + pp.beta);E= spdiags(e,0,G,G);

      % Define derivatives
      dFc_ds = (speye(G,G) - X*inv(E)*F)*L;
      dFc_dx = -E;
      dFf_ds = J'*J;
      dFf_dx = L';

      dsdx = -[dFc_ds, dFc_dx; dFf_ds, dFf_dx] \ ...
              [ f-E*x; J'*(J*s-d) + L'*x ];

      ds = dsdx(1:N);
      dx = x_update(x, dsdx(N+(1:G)));

      s= s + ds; x= x + dx;
fprintf('+');
   end

function s= pdipm_1_1( J,W,L,d, pp);
   [M,N] = size(J); % M measurements, N parameters
   [D  ] = size(L,1); % E edges
   y= zeros( D, 1 ); % dual var - start with zeros
   s= zeros( N, 1 ); % solution - start with zeros
   x= zeros( M, 1 ); % dual var - start with zeros

   for loop = 1:pp.max_iter
      % Define variables
      g = L*s;                 G= spdiags(g,0,D,D);
      r = sqrt(g.^2 + pp.beta);R= spdiags(r,0,D,D); % S in paper
                               Y= spdiags(y,0,D,D);

      f = J*s - d;             F= spdiags(f,0,M,M);
      e = sqrt(f.^2 + pp.beta);E= spdiags(e,0,M,M);
                               X= spdiags(x,0,M,M);

      % Define derivatives
      As1 = sparse(N,N);
      As2 = (speye(M,M) - X*inv(E)*F) * J;
      As3 = (speye(D,D) - Y*inv(R)*G) * L;
      Ax1 = J'*W;
      Ax2 = -E;
      Ax3 = sparse(D,M);
      Ay1 = L';
      Ay2 = sparse(M,D);
      Ay3 = -R;
      B1  = J'*W*x + L'*y;
      B2  = f - E*x;
      B3  = g - R*y;

      DD = -[As1,Ax1,Ay1; ...
             As2,Ax2,Ay2; ...
             As3,Ax3,Ay3] \ [B1;B2;B3];

      ds = DD(1:N);
      dx = x_update(x, DD(N+(1:M)));
      dy = x_update(y, DD(N+M+(1:D)));

      s= s + ds;
      x= x + dx;
      y= y + dy;
fprintf('+');
   end

% abs(x + dx) must be <= 1
function dx = x_update( x, dx)
   dx(dx==0) = eps; % can't have zeros
   sx = sign(dx);
   % space to limits in direction of x
   clr = sx - x;
   % how much to multiply by to get to limits
   fac = clr./dx;
   % choose min amount to get to limits
   dx = dx*min(fac);

function pp= process_parameters(imdl);
   try    pp.max_iter = imdl.parameters.max_iterations;
   catch  pp.max_iter = 10;
   end

   try    pp.min_change = imdl.parameters.min_change;
   catch  pp.min_change = 0;
   end

   try    pp.beta = imdl.ab_pdipm.beta; 
   catch  pp.beta = 1e-8;
   end

   try    pp.norm_data = imdl.ab_pdipm.norm_data;
   catch  pp.norm_data = 2;
   end

   try    pp.norm_image = imdl.ab_pdipm.norm_image;
   catch  pp.norm_image = 2;
   end


