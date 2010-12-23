function img=pdipm_abs( inv_model, data);
% PDIPM_ABS  inverse solver for absolute data using Primal/Dual interior point method
% img= pdipm_abs( inv_model, data);
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data       => vector of eit data
%
%  inv_model.pdipm_abs.norm_data  1 or 2 (DEFAULT 2)
%  inv_model.pdipm_abs.norm_prior 1 or 2 (DEFAULT 2)
%  inv_model.pdipm_abs.beta     (default 1e-6)
%
% Parameters:
%  max_iter =  inv_model.parameters.max_iteration (default 10)
%      Max number of iterations before stopping
%  min change = inv_model.parameters.min_change   (default 0)
%      Min Change in objective fcn (norm(y-Jx)^2 + hp*TV(x)) before stopping
% beta is the parameter that smooths the TV functional

% (C) 2010 Andrea Borsic + Andy Adler. License: GPL v2 or v3
% $Id$


pp= process_parameters(inv_model);

img_bkgnd = homogeneous_estimate( inv_model, data );

alpha=calc_hyperparameter( inv_model );
L= calc_R_prior( inv_model );
W= calc_meas_icov( inv_model );

img= feval(pp.fn, img_bkgnd, W,alpha*L,data, pp);

img.name = sprintf('pdipm_abs-nd%d-ni%d',pp.norm_data,pp.norm_image);

% This is the Gauss-Newton algorithm
%   for the linear case it is: s= (J'*W*J + L'*L)\J'*W*d;
function img= pdipm_2_2(  img,W,L,d, pp);
   img0 = img;
   hp2RtR = L'*L;
   for i = 1:pp.max_iter
     dv = sim_diff( img, d);
     J = calc_jacobian( img );

     RDs = hp2RtR*(img0.elem_data - img.elem_data);
     ds = (J'*W*J + hp2RtR)\(J'*dv + RDs);

     img = line_optimize(img, ds, d);

     loop_display(i)
   end

function img= pdipm_1_2( img,W,L,d, pp);
   [M]   = size(W,1); % M measurements
   [jnk,N] = size(L); % E edges, N parameters
   x= zeros( M, 1 ); % dual var - start with zeros

   for loop = 1:pp.max_iter
      % Jacobian
      J = calc_jacobian( img );
      dv = -sim_diff(img, d);
      s = img.elem_data;

      % Define variables
      f = dv;                  F= spdiags(f,0,M,M);
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
      img = line_optimize(img, ds, d);

      dx = x_update(x, dsdx(N+(1:M)));
      x= x + dx;

      loop_display(i)
   end

function img= pdipm_2_1(img,W,L,d, pp);
   [M]   = size(W,1); % M measurements
   [G,N] = size(L); % E edges, N parameters
   x= zeros( G, 1 ); % dual var - start with zeros

   for loop = 1:pp.max_iter
      % Jacobian
      J = calc_jacobian( img );
      dv = -sim_diff(img, d);
      
      % Define variables
      f = L*img.elem_data;     F= spdiags(f,0,G,G);
                               X= spdiags(x,0,G,G);
      e = sqrt(f.^2 + pp.beta);E= spdiags(e,0,G,G);

      % Define derivatives
      dFc_ds = (speye(G,G) - X*inv(E)*F)*L;
      dFc_dx = -E;
      dFf_ds = J'*J;
      dFf_dx = L';

      dsdx = -[dFc_ds, dFc_dx; dFf_ds, dFf_dx] \ ...
              [ f-E*x; J'*dv + L'*x ];

      ds = dsdx(1:N);
      img = line_optimize(img, ds, d);

      dx = x_update(x, dsdx(N+(1:G)));
      x= x + dx;

      loop_display(i)
   end

%   img0 = img;
%   hp2RtR = L'*L;
%   for i = 1:pp.max_iter
%     vsim = fwd_solve( img ); 
%     dv = calc_difference_data( vsim , d, img.fwd_model);
%     J = calc_jacobian( img );
%
%     RDx = hp2RtR*(img0.elem_data - img.elem_data);
%     dx = (J'*W*J + hp2RtR)\(J'*dv + RDx);
%
%     img = line_optimize(img, dx, d);
%
%     loop_display(i)
%   end

function img= pdipm_1_1( img,W,L,d, pp);
   [M]   = size(W,1); % M measurements
   [D,N] = size(L); % E edges, N parameters
   s= img.elem_data;
   y= zeros( D, 1 ); % dual var - start with zeros
   x= zeros( M, 1 ); % dual var - start with zeros

   for loop = 1:pp.max_iter
      % Jacobian
      J = calc_jacobian( img );
      dv = -sim_diff(img, d);

      % Define variables
      g = L*img.elem_data;     G= spdiags(g,0,D,D);
      r = sqrt(g.^2 + pp.beta);R= spdiags(r,0,D,D); % S in paper
                               Y= spdiags(y,0,D,D);

      f = dv;                  F= spdiags(f,0,M,M);
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
      img = line_optimize(img, ds, d);

      dx = x_update(x, DD(N+(1:M)));
      x= x + dx;

      dy = x_update(y, DD(N+M+(1:D)));
      y= y + dy;

      loop_display(i)
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

   try    pp.beta = imdl.pdipm_abs.beta; 
   catch  pp.beta = 1e-8;
   end

   try    pp.norm_data = imdl.pdipm_abs.norm_data;
   catch  pp.norm_data = 2;
   end

   try    pp.norm_image = imdl.pdipm_abs.norm_image;
   catch  pp.norm_image = 2;
   end

   if     pp.norm_data==2 && pp.norm_image==2;
      pp.fn = @pdipm_2_2;
   elseif pp.norm_data==2 && pp.norm_image==1;
      pp.fn = @pdipm_2_1;
   elseif pp.norm_data==1 && pp.norm_image==2;
      pp.fn = @pdipm_1_2;
   elseif pp.norm_data==1 && pp.norm_image==1;
      pp.fn = @pdipm_1_1;
   else
      error('norm_data and norm_image should be 1 or 2');
   end



% Fit a parabola to the linefit and pick the best point
% This is faster than doing an exhaustive search
function  img = line_optimize(imgk, dx, data1);
  flist = [ 0.1,  0.5, 1.0];
  clim = mean(imgk.elem_data)/10; % prevent zero and negative conductivity
  img = imgk;
  for i = 1:length(flist);
     img.elem_data = imgk.elem_data + flist(i)*dx;
     img.elem_data(img.elem_data <= clim ) = clim;
     dv = sim_diff( img, data1);
     mlist(i) = norm(dv);
  end
  pf = polyfit(flist, mlist, 2);
  fmin = -pf(2)/pf(1)/2; % poly minimum
  fmin(fmin>1) = 1; fmin(fmin<0) = 0;

  img.elem_data = imgk.elem_data + flist(i)*dx;
  img.elem_data(img.elem_data <= clim ) = clim;

function img = homogeneous_estimate( imdl, data );
   img = calc_jacobian_bkgnd( imdl );
   vs = fwd_solve(img);
   data = data_vector( data, imdl );

   pf = polyfit(data,vs.meas,1);

   img.elem_data = img.elem_data*pf(1);

function data = data_vector( data, imdl );
   if isstruct(data)
      data = data.meas;
   else
     meas_select = [];
     try
        meas_select = imdl.fwd_model.meas_select;
     end
     if length(data) == length(meas_select)
        data = data(meas_select);
     end
   end

function dv = sim_diff( img, data1);
  vsim = fwd_solve( img );
  dv = calc_difference_data( vsim , data1, img.fwd_model);

function loop_display(i)
   fprintf('+');
