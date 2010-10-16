function img=pdipm_diff( inv_model, data1, data2)
% PDIPM_DIFF inverse solver for difference data using Primal/Dual interior point method
% img= ab_pdipm( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
%  inv_model.pdipm_diff.norm_data  1 or 2 (DEFAULT 2)
%  inv_model.pdipm_diff.norm_image 1 or 2 (DEFAULT 2)
%  inv_model.pdipm_diff.beta     (default 1e-6)
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
if pp.norm_data==1
  W = sqrt(W); % sW is in units of volts
end

if     pp.norm_data==2 && pp.norm_image==2
  x= pdipm_2_2( J,W,alpha*L,d, pp);
elseif pp.norm_data==2 && pp.norm_image==1
  x= pdipm_2_1( J,W,alpha*L,d, pp);
elseif pp.norm_data==1 && pp.norm_image==2
  x= pdipm_1_2( J,W,alpha*L,d, pp);
elseif pp.norm_data==1 && pp.norm_image==1
  x= pdipm_1_1( J,W,alpha*L,d, pp);
end

% create a data structure to return
img.name = 'pdipm_diff';
img.elem_data = x;
img.fwd_model = fwd_model;

function s= pdipm_2_2( J,W,L,d, pp);
   [s]= initial_values( J, L, pp);

   R = L'*L;
   ds= (J'*W*J + R)\(J'*W*(d - J*s) - R*s);
   s= s + ds;

function s= pdipm_1_2( J,W,L,d, pp);
   [s,x,jnk,sz]= initial_values( J, L, pp);

   I_M = speye(sz.M, sz.M);
   for loop = 1:pp.max_iter
      % Define variables
      f = J*s - d;             F= spdiag(f);
                               X= spdiag(x);
      e = sqrt(f.^2 + pp.beta);E= spdiag(e);

      % Define derivatives
      dFc_ds = (I_M - X*inv(E)*F)*J;
      dFc_dx = -E;
      dFf_ds = L'*L;
      dFf_dx = J'*W;

      dsdx = -[dFc_ds, dFc_dx; dFf_ds, dFf_dx] \ ...
              [ f-E*x; J'*W*x + L'*L*s ];

      ds =             dsdx(      1:sz.N);
      dx = x_update(x, dsdx(sz.N+(1:sz.M)));

      s= s + ds; x= x + dx;
fprintf('+');
   end

function s= pdipm_2_1( J,W,L,d, pp);
   [s,jnk,x,sz]= initial_values( J, L, pp);

   I_D = speye(sz.D, sz.D);
   for loop = 1:pp.max_iter
      % Define variables
      f = L*s;                 F= spdiag(f);
                               X= spdiag(x);
      e = sqrt(f.^2 + pp.beta);E= spdiag(e);

      % Define derivatives
      dFc_ds = (I_D - X*inv(E)*F)*L;
      dFc_dx = -E;
%     dFf_ds = J'*J;
      dFf_ds = J'*W*J;
      dFf_dx = L';

      dsdx = -[dFc_ds, dFc_dx; dFf_ds, dFf_dx] \ ...
              [ f-E*x; J'*W*(J*s-d) + L'*x ];
%             [ f-E*x; J'*(J*s-d) + L'*x ];

      ds =             dsdx(      1:sz.N );
      dx = x_update(x, dsdx(sz.N+(1:sz.D)));

      s= s + ds; x= x + dx;
fprintf('+');
   end

function s= pdipm_1_1( J,W,L,d, pp);
   [s,x,y,sz]= initial_values( J, L, pp);

   I_M = speye(sz.M,sz.M); 
   I_D = speye(sz.D,sz.D); 
   Z_N = sparse(sz.N,sz.N);
   Z_DM= sparse(sz.D,sz.M);
   for loop = 1:pp.max_iter
      % Define variables
      g = L*s;                 G= spdiag(g);
      r = sqrt(g.^2 + pp.beta);R= spdiag(r); % S in paper
                               Y= spdiag(y);

      f = J*s - d;             F= spdiag(f);
      e = sqrt(f.^2 + pp.beta);E= spdiag(e);
                               X= spdiag(x);

      % Define derivatives
      As1 = Z_N;
      As2 = (I_M - X*inv(E)*F) * J;
      As3 = (I_D - Y*inv(R)*G) * L;
      Ax1 = J'*W;
      Ax2 = -E;
      Ax3 = Z_DM;
      Ay1 = L';
      Ay2 = Z_DM';
      Ay3 = -R;
      B1  = J'*W*x + L'*y;
      B2  = f - E*x;
      B3  = g - R*y;

      DD = -[As1,Ax1,Ay1; ...
             As2,Ax2,Ay2; ...
             As3,Ax3,Ay3] \ [B1;B2;B3];

      ds = DD(1:sz.N);
      dx = x_update(x, DD(sz.N +        (1:sz.M)) );
      dy = x_update(y, DD(sz.N + sz.M + (1:sz.D)) );

      s= s + ds;
      x= x + dx;
      y= y + dy;
fprintf('+');
   end

% fix matlab's stupid verbose spdiags function
function sM = spdiag(V)
   lV = length(V);
   sM = spdiags(V,0,lV,lV);

function [s,x,y,sz]= initial_values( J, L, pp);
   [sz.M,sz.N] = size(J); % M measurements, N parameters
   [sz.D  ] = size(L,1); % E edges
   y= zeros( sz.D, 1 ); % dual var - start with zeros
   s= zeros( sz.N, 1 ); % solution - start with zeros
   x= zeros( sz.M, 1 ); % dual var - start with zeros

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

   try    pp.beta = imdl.pdipm_diff.beta; 
   catch  pp.beta = 1e-8;
   end

   try    pp.norm_data = imdl.pdipm_diff.norm_data;
   catch  pp.norm_data = 2;
   end

   try    pp.norm_image = imdl.pdipm_diff.norm_image;
   catch  pp.norm_image = 2;
   end
