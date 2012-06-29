function img= GN_absLogCond_solve( inv_model, data1);
% INV_SOLVE_ABS_GNR absolute solver using Gauss Newton approximation
% img= gn_abs_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => EIT data
%
% Parameters:
%   inv_model.parameters.max_iterations = N_max iter

% (C) 2012 Nolwenn Lesparre. License: GPL version 2 or version 3
% $Id$

% Step 1: fit to background
img = homogeneous_estimate( inv_model, data1 );
hp  = calc_hyperparameter( inv_model );
RtR = calc_RtR_prior( inv_model );
% RtR= speye(size(img.elem_data,1));
W   = calc_meas_icov( inv_model );
hp2RtR= hp*RtR;

iters = 1;
try 
   iters = inv_model.parameters.max_iterations; 
end
img0 = img;
img0.logRes= -log(img0.elem_data);
for i = 1:iters  
  vsim = fwd_solve( img ); 
  dv = calc_difference_data( vsim , data1, img.fwd_model);
  dv= inv_model.parameters.normalisation*dv;
  img.logRes= -log(img.elem_data);
  dCond_dlogRes= -exp(-img.logRes);
  
  RDx = hp2RtR*(img0.logRes - img.logRes);
  J = calc_jacobian( img );
  J = J.*repmat((dCond_dlogRes),1,size(data1,1))';
  J= inv_model.parameters.normalisation*J;
  dx = (J'*W*J + hp2RtR)\(J'*W*dv + RDx);
 
  img = line_optimize(img, dx, data1);
end


% Fit a parabola to the linefit and pick the best point
% This is better than doing an exhaustive search
function  img = line_optimize(imgk, dx, data1);
    img = imgk;
%    flist =  logspace(-7,0,15);%0.01:0.01:0.15;
%    mlist= flist*0;
% % %   clim = mean(imgk.elem_data)/10; % prevent zero and negative conductivity
%   for i = 1:length(flist);
% %      img.elem_data = imgk.elem_data + flist(i)*dx;
% %      img.elem_data(img.elem_data <= clim ) = clim;
%      img.logRes= imgk.logRes + flist(i)*dx;
%      img.elem_data= exp(-img.logRes);
%      vsim = fwd_solve( img );
%      dv = calc_difference_data( vsim , data1, img.fwd_model);
%      mlist(i) = norm(dv);
%   end
% %   pf = polyfit(flist, mlist, 2);
% %   fmin = -pf(2)/pf(1)/2; % poly minimum
%   [minm,iminm]= min(mlist);
%   fmin= flist(iminm);
  fmin= 0.01;
%   figure; loglog(flist,mlist,'xr',flist*0+fmin,mlist)
%   fmin(fmin>1) = 1; fmin(fmin<0) = 0;  
%   fmin= 0.05;
  img.logRes= imgk.logRes + fmin*dx;
  img.elem_data= exp(-img.logRes);
  
%   img.elem_data = imgk.elem_data + fmin*dx;
%   img.elem_data(img.elem_data <= clim ) = clim;
  img.res_data= exp(img.logRes);
  
  
function img = homogeneous_estimate( imdl, data );
   img = calc_jacobian_bkgnd( imdl );
   vs = fwd_solve(img);

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

   pf = polyfit(data,vs.meas,1);

   if isfield(img.fwd_model,'coarse2fine');
% TODO: the whole coarse2fine needs work here.
%   what happens if c2f doesn't cover the whole region
      nc = size(img.fwd_model.coarse2fine,2);
      img.elem_data = mean(img.elem_data)*ones(nc,1)*pf(1);
   else
      img.elem_data = img.elem_data*pf(1);
   end
