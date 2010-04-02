function img= GN_abs_solve( inv_model, data1);
% GN_ABS_SOLVER absolute solver using Gauss Newton approximation
% img= gn_abs_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => EIT data
%
% Parameters:
%   inv_model.parameters.max_iteration = N_max iter

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

% Step 1: fit to background
img = homogeneous_estimate( inv_model, data1 );

hp  = calc_hyperparameter( inv_model );
RtR = calc_RtR_prior( inv_model );
W   = calc_meas_icov( inv_model );
hp2RtR= hp*RtR;

iters = 1;
try 
   iters = inv_model.parameters.max_iteration;
end

img0 = img;
for i = 1:iters; 
  vsim = fwd_solve( img ); 
  dv = calc_difference_data( vsim , data1, img.fwd_model);
  J = calc_jacobian( img );

  RDx = hp2RtR*(img0.elem_data - img.elem_data);
  dx = (J'*W*J + hp2RtR)\(J'*dv - RDx);

  img = line_optimize(img, dx, data1);
end

% Fit a parabola to the linefit and pick the best point
% This is better than doing an exhaustive search
function  img = line_optimize(imgk, dx, data1);
  flist = [ 0.1,  0.5, 0.9];
  clim = mean(imgk.elem_data)/10; % prevent zero and negative conductivity
  img = imgk;
  for i = 1:length(flist);
     img.elem_data = imgk.elem_data + flist(i)*dx;
     img.elem_data(img.elem_data <= clim ) = clim;
     vsim = fwd_solve( img );
     dv = calc_difference_data( vsim , data1, img.fwd_model);
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

   img.elem_data = img.elem_data*pf(1);
