function img= inv_solve_abs_GN( inv_model, data1);
% INV_SOLVE_ABS_GNR absolute solver using Gauss Newton approximation
% img= gn_abs_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => EIT data
%
% Parameters:
%   inv_model.parameters.max_iterations = N_max iter

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
   iters = inv_model.parameters.max_iterations; 
end

img0 = physics_data_mapper(img);
if ~strcmp(img0.current_physics,'conductivity')
   error('inv_solve_abs_GN does not work for %s',img0.current_physics);
end
for i = 1:iters  
  vsim = fwd_solve( img ); 
  dv = calc_difference_data( vsim , data1, img.fwd_model);
  J = calc_jacobian( img );

  % create elem_data
  tmp = physics_data_mapper(img);
  RDx = hp2RtR*(img0.elem_data - tmp.elem_data);
  dx = (J'*W*J + hp2RtR)\(J'*dv + RDx);

  img = line_optimize(img, dx, data1);
end

% Fit a parabola to the linefit and pick the best point
% This is better than doing an exhaustive search
function  img = line_optimize(imgk, dx, data1)
  % get elem_data
  imgk = physics_data_mapper(imgk);
  flist = [ 0.1,  0.5, 1.0];
  clim = mean(imgk.elem_data)/10; % prevent zero and negative conductivity
  img = imgk;
  for i = 1:length(flist);
     img.elem_data = imgk.elem_data + flist(i)*dx;
     img.elem_data(img.elem_data <= clim ) = clim;
     vsim = fwd_solve( physics_data_mapper( img, 1 ) );
     dv = calc_difference_data( vsim , data1, img.fwd_model);
     mlist(i) = norm(dv);
  end
  pf = polyfit(flist, mlist, 2);
  fmin = -pf(2)/pf(1)/2; % poly minimum for a 2nd order poly
  fmin(fmin>1) = 1; fmin(fmin<0) = 0; % set to limits of [0,1]

  img.elem_data = imgk.elem_data + fmin*dx;
  img.elem_data(img.elem_data <= clim ) = clim;
  
  % recreate physics
  img = physics_data_mapper(img,1);
  

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

   % create elem_data
   img = physics_data_mapper(img);
   
   if isfield(img.fwd_model,'coarse2fine');
% TODO: the whole coarse2fine needs work here.
%   what happens if c2f doesn't cover the whole region
      nc = size(img.fwd_model.coarse2fine,2);
      img.elem_data = mean(img.elem_data)*ones(nc,1)*pf(1);
   else
      img.elem_data = img.elem_data*pf(1);
   end
   
   % remove elem_data
   img = physics_data_mapper(img,1);
