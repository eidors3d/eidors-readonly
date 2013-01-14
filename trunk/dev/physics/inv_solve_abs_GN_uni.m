function img= inv_solve_abs_GN( inv_model, data1);
% INV_SOLVE_ABS_GNR absolute solver using Gauss Newton approximation
% img= gn_abs_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => EIT data
%
% Parameters:
%   inv_model.parameters.max_iterations = N_max iter
% Parameters: (will override parameters field)
%   Maximum Iterations
%    inv_model.inv_solve_abs_GN.max_iterations = N_max iter
%   Line Optimize parameters
%    inv_model.inv_solve_abs_GN.line_optimize.perturb
%    inv_model.inv_solve_abs_GN.line_optimize.min_value
%    inv_model.inv_solve_abs_GN.line_optimize.max_value
%   Perform initial estimate of Homogeneous background
%    inv_model.inv_solve_abs_GN.do_starting_estimate  (default 1)

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

% Step 1: fit to background
img = homogeneous_estimate( inv_model, data1 );

hp  = calc_hyperparameter( inv_model );
RtR = calc_RtR_prior( inv_model );
W   = calc_meas_icov( inv_model );
hp2RtR= hp*RtR;

opt = parse_options(inv_model);
iters = opt.max_iter;

img0 = physics_data_mapper(img);

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


function img = homogeneous_estimate( imdl, data )
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

% TODO: the two cases are very different. c2f case should match other
      nc = size(img.fwd_model.coarse2fine,2);
      img.elem_data = mean(img.elem_data)*ones(nc,1)*pf(1);
   else
      img.elem_data = img.elem_data*pf(1);
   end
   
   % remove elem_data
   img = physics_data_mapper(img,1);
   
   
   
function opt = parse_options(imdl)

opt.max_iter = 1;
try
   opt.max_iter = imdl.parameters.max_iterations;
end
if isfield(imdl, 'inv_solve_abs_GN');
   opt = imdl.inv_solve_abs_GN;
end
