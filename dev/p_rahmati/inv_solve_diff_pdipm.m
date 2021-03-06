function img=inv_solve_diff_pdipm( inv_model, data1, data2,varargin)
% INV_SOLVE_DIFF_PDIPM inverse solver for difference data using Primal/Dual interior point method
% img= ab_pdipm( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
%  inv_model.inv_solve_diff_pdipm.norm_data  1 or 2 (DEFAULT 2)
%  inv_model.inv_solve_diff_pdipm.norm_image 1 or 2 (DEFAULT 2)
%  inv_model.inv_solve_diff_pdipm.beta     (default 1e-6)
%
% Parameters:
%  max_iters =  inv_model.parameters.max_iterations (default 10)
%      Max number of iterations before stopping
%  min change = inv_model.parameters.min_change   (default 0)
%      Min Change in objective fcn (norm(y-Jx)^2 + hp*TV(x)) before stopping
% beta is the parameter that smooths the TV functional

% (C) 2008 Andrea Borsic. License: GPL version 2 or version 3
% $Id$

%-----------------------------Peyman's code -----------------
%if pp.norm_data_1==1 && pp.norm_image_1==1 && pp.norm_data_2==2 && pp.norm_image_2==2%L1L1L2L2 
switch nargin%('inv_solve_diff_pdipm')
        case 3
           zeta=.6;%default
           eta=.6; % default
           fprintf('zeta= eta=0.6 is selected as the default values for the L1L1L2L2!');
   
    otherwise 
            zeta= varargin{1};
            eta=varargin{2};
            fprintf('zeta and eta is selected by user for the L1L1L2L2!');

end

%end
%--------------------------------------------------------------------------------------
pp= process_parameters(inv_model);

fwd_model= inv_model.fwd_model;

d=calc_difference_data( data1, data2, fwd_model);

img_bkgnd=calc_jacobian_bkgnd( inv_model );
%J=calc_jacobian( fwd_model, img_bkgnd);
%J=calc_move_jacobian( fwd_model, img_bkgnd);
alpha=calc_hyperparameter( inv_model );
%-------calc RtR prior or matrix L---------
movement_compensation=0;
%L=calc_R_prior( inv_model );
if movement_compensation
   L = laplace_image_prior(inv_model);
% L = gaussian_HPF_prior(inv_model);
%    idx= 901:1024; % Outer row
%    idx= 785:1024; % Two outer rows
%    RtR(idx,idx) = 4*RtR(idx,idx);
 %HP = 0.06;
       P2= 1e-2; % NT = 0.864
%      HP = 0.10; P2= 1e-1; % NF = 0.506
       inv_model.fwd_model.normalize_measurements=1;
       inv_model.fwd_model.jacobian = @calc_move_jacobian;
       inv_model.RtR_prior =          @aa_e_move_image_prior;
       inv_model.aa_e_move_image_prior.parameters = sqrt(P2/1); 
%       imdl.hyperparameter.value = HP;
       inv_model.inv_solve.select_parameters = 1:size(inv_model.fwd_model.elems,1);

       inv_model.aa_e_move_image_prior.RegC.func = @gaussian_HPF_prior;
       inv_model.aa_e_move_image_prior.RegC.func = @d2_outer_weight_laplace_prior;
%----------------------------------
     %J=calc_jacobian( fwd_model, img_bkgnd);
 J=jacobian_movement(fwd_model, img_bkgnd);
 J=J(:,1:size(fwd_model.elems,1));

else 
    J=calc_jacobian( fwd_model, img_bkgnd);
    L=calc_R_prior( inv_model );
end
W= calc_meas_icov( inv_model );
%--------------Peyman's code-----------
if pp.norm_data_1==1
  W = sqrt(W); % sW is in units of volts
end
%----------------------- Generalized PDIPM parameters (Peyman's code)------
% Norms           |   L2L2 | L1L2 | L2L1 | L1L1 | L1L1L2L2
% pp.norm_data_1  |    0   |   1  | 0    | 1    | 1
% pp.norm_image_1 |    0   |   0  | 1    | 1    | 1
% pp.norm_data_2  |    2   |   0  | 2    | 0    | 2
% pp.norm_image_2 |    2   |   2  | 0    | 0    | 2
% \zeta           |    0   |   1  | 0    | 1    | 0<\zeta<1
% \eta            |    0   |   0  | 1    | 1    | 0<\eta<1
%----------------------------------------------------------
varargin
if     pp.norm_data_2==2 && pp.norm_image_2==2 && pp.norm_data_1==0 && pp.norm_image_1==0 %L2L2
  x= pdipm_2_2( J,W,alpha*L,d, pp);
elseif pp.norm_data_2==2 && pp.norm_image_1==1&& pp.norm_data_1==0 && pp.norm_image_2==0%L2L1
  x= pdipm_2_1( J,W,alpha*L,d, pp);
elseif pp.norm_data_1==1 && pp.norm_image_2==2 && pp.norm_data_2==0 && pp.norm_image_1==0%L1L2
  x= pdipm_1_2( J,W,alpha*L,d, pp);
elseif pp.norm_data_1==1 && pp.norm_image_1==1 && pp.norm_data_2==0 && pp.norm_image_2==0%L1L1
  x= pdipm_1_1( J,W,alpha*L,d, pp);
elseif pp.norm_data_1==1 && pp.norm_image_1==1 && pp.norm_data_2==2 && pp.norm_image_2==2%L1L1L2L2
  x= pdipm_1_1_2_2( J,W,alpha*L,d, pp,zeta, eta);
end

% create a data structure to return
img.name = 'inv_solve_diff_pdipm';
img.type='image';
img.elem_data = x;
img.fwd_model = fwd_model;
function RtR = d2_outer_weight_laplace_prior( inv_model);
   RtR = laplace_image_prior(inv_model);
function s= pdipm_2_2( J,W,L,d, pp);
   [s]= initial_values( J, L, pp);

   R = L'*L;
   ds= (J'*W*J + R)\(J'*W*(d - J*s) - R*s);
   s= s + ds;

function m= pdipm_1_2( J,W,L,d, pp);
   [m,x,jnk,sz]= initial_values( J, L, pp);

   I_M = speye(sz.M, sz.M);
   for loop = 1:pp.max_iter
      % Define variables
      f = J*m - d;             F= spdiag(f);
                               X= spdiag(x);
      e = sqrt(f.^2 + pp.beta);E= spdiag(e);

      % Define derivatives
      dFc_dm = (I_M - X*inv(E)*F)*J;
      dFc_dx = -E;
      dFf_dm = L'*L; % SHOULD NOT THIS BE 2*L'*L ?! 
      dFf_dx = J'*W;

      dmdx = -[dFc_dm, dFc_dx; dFf_dm, dFf_dx] \ ...
              [ f-E*x; J'*W*x + L'*L*m ];% SHOULD NOT THIS BE 2*L'*L ?!

      dm =             dmdx(      1:sz.N);
      dx = x_update(x, dmdx(sz.N+(1:sz.M)));

      m= m + dm; x= x + dx;
      loop_display(i)
debug([mean(abs([m,dm])) mean(abs([x,dx]))])
      pp = manage_beta(pp);
   end

function m= pdipm_2_1( J,W,L,d, pp);
   [m,jnk,y,sz]= initial_values( J, L, pp);

   I_D = speye(sz.D, sz.D);
   for loop = 1:pp.max_iter
      % Define variables
      g = L*m;                 G= spdiag(g);
                               Y= spdiag(y);
      s = sqrt(g.^2 + pp.beta);S= spdiag(s);

      % Define derivatives
      dFf_dm = 2*J'*W*J;% SHOULD NOT THIS BE 2*J'*(W'*W)*J?!
      dFf_dy = L';
      dFc_dm = (I_D - Y*inv(S)*G)*L;
      dFc_dy = -S;

      dmdy = -[dFf_dm, dFf_dy; dFc_dm, dFc_dy] \ ...
              [ J'*W*(J*m-d) + L'*y; g-S*y ];% SHOULD NOT THIS BE 2*J'*W*(J*m-d) + L'*y?!

      dm =             dmdy(      1:sz.N );
      dy = x_update(y, dmdy(sz.N+(1:sz.D)));

      m= m + dm; y= y + dy;
      loop_display(i)
debug([mean(abs([m,dm])), mean(abs([y,dy]))]);
      pp = manage_beta(pp);
   end

function m= pdipm_1_1( J,W,L,d, pp);
   [m,x,y,sz]= initial_values( J, L, pp);

   I_M = speye(sz.M,sz.M); 
   I_D = speye(sz.D,sz.D); 
   Z_N = sparse(sz.N,sz.N);
   Z_DM= sparse(sz.D,sz.M);
   for loop = 1:pp.max_iter
      % Define variables
      g = L*m;                 G= spdiag(g);
      r = sqrt(g.^2 + pp.beta);R= spdiag(r); % S in paper
                               Y= spdiag(y);

      f = J*m - d;             F= spdiag(f);
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

      dm = DD(1:sz.N);
      dx = x_update(x, DD(sz.N +        (1:sz.M)) );
      dy = x_update(y, DD(sz.N + sz.M + (1:sz.D)) );

      m= m + dm;
      x= x + dx;
      y= y + dy;
      loop_display(i)
debug([mean(abs([m,dm])), mean(abs([x,dx])), mean(abs([y,dy]))]);
      pp = manage_beta(pp);
   end
function m= pdipm_1_1_2_2( J,W,L,d, pp, zeta, eta);
   [m,x,y,sz]= initial_values( J, L, pp);

   I_M = speye(sz.M,sz.M); 
   I_D = speye(sz.D,sz.D); 
   Z_N = sparse(sz.N,sz.N);
   Z_DM= sparse(sz.D,sz.M);
   for loop = 1:pp.max_iter
      % Define variables
      g = L*m;                 G= spdiag(g);
      r = sqrt(g.^2 + pp.beta);R= spdiag(r); % S in paper
                               Y= spdiag(y);

      f = J*m - d;             F= spdiag(f);
      e = sqrt(f.^2 + pp.beta);E= spdiag(e);
                               X= spdiag(x);

      % Define derivatives
      As1 = 2*(1-zeta)*J'*W*J+2*(1-eta)*(L'*L);%Z_N; or 2*(1-zeta)*J'*(W'*W)*J+2*(1-eta)*(L'*L);
      As2 = (I_M - X*inv(E)*F) * J;
      As3 = (I_D - Y*inv(R)*G) * L;
      Ax1 = zeta*J'*W;
      Ax2 = -E;
      Ax3 = Z_DM;
      Ay1 = eta*L';
      Ay2 = Z_DM';
      Ay3 = -R;
      B1  = zeta*J'*W*x + eta* L'*y+2*(1-zeta)*J'*W*f+2*(1-eta)*L'*g;%zeta*J'*W*x + eta* L'*y+2*(1-zeta)*J'*(W'*W)*F+2*(1-eta)*L'*G;
      B2  = f - E*x;
      B3  = g - R*y;

      DD = -[As1,Ax1,Ay1; ...
             As2,Ax2,Ay2; ...
             As3,Ax3,Ay3] \ [B1;B2;B3];

      dm = DD(1:sz.N);
      dx = x_update(x, DD(sz.N +        (1:sz.M)) );
      dy = x_update(y, DD(sz.N + sz.M + (1:sz.D)) );

      m= m + dm;
      x= x + dx;
      y= y + dy;
      loop_display(i)
debug([mean(abs([m,dm])), mean(abs([x,dx])), mean(abs([y,dy]))]);
      pp = manage_beta(pp);
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
   %s= zeros( size(L,2), 1 ); % solution - start with zeros
   x= zeros( sz.M, 1 ); % dual var - start with zeros

% abs(x + dx) must be <= 1
function dx = x_update( x, dx)
 % can't have zeros
   dx(dx==0) = eps;
 % space to limits in direction of dx
   sx = sign(dx);
   clr = sx - x;
   % how much to multiply by to get to limits
   fac = clr./dx;
   % choose min amount to get to limits
   dx = dx*min(abs(fac));
%  dx = dx*0.9;

function debug(vals)
%  disp(vals)

function pp = manage_beta(pp);
   pp.beta = pp.beta * pp.beta_reduce;
   if pp.beta < pp.beta_minimum;
      pp.beta = pp.beta_minimum;
   end

function pp= process_parameters(imdl);
   try    pp.max_iter = imdl.parameters.max_iterations;
   catch  pp.max_iter = 10;
   end

   try    pp.min_change = imdl.parameters.min_change;
   catch  pp.min_change = 0;
   end

   try    pp.beta = imdl.inv_solve_diff_pdipm.beta; 
   catch  pp.beta = 1e-6;
   end

   pp.beta_reduce = 0.2;
   pp.beta_minimum= 1e-16;
%-----------------------peyman's code (default is L2L2)--------------------
   try    pp.norm_data_2 = imdl.inv_solve_diff_pdipm.norm_data_2;
   catch  pp.norm_data_2 = 2;
   end

   try    pp.norm_image_2 = imdl.inv_solve_diff_pdipm.norm_image_2;
   catch  pp.norm_image_2 = 2;
   end
   
   try    pp.norm_data_1 = imdl.inv_solve_diff_pdipm.norm_data_1;
   catch  pp.norm_data_1 = 0;
   end

   try    pp.norm_image_1 = imdl.inv_solve_diff_pdipm.norm_image_1;
   catch  pp.norm_image_1 = 0;
   end
   %-----------------------------------------------------

function loop_display(i)
   fprintf('+');
