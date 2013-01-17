function img= inv_solve_abs_GN( inv_model, data0);
%INV_SOLVE_ABS_GN Absolute solver using Gauss Newton approximation
% img= inv_solve_abs_GN( inv_model, data0)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data0      => EIT data
%
% Parameters:
%   inv_model.parameters.max_iterations = N_max iter
% Parameters: (will override parameters field)
%   Maximum Iterations
%    inv_model.inv_solve_abs_GN.max_iterations = N_max iter
%   Line Optimize parameters
%    inv_model.inv_solve_abs_GN.line_optimize.perturb
%   Limits on solution
%    inv_model.inv_solve_abs_GN.min_value
%    inv_model.inv_solve_abs_GN.max_value
%   Update function after each iteration (more details below)
%    inv_model.inv_solve_abs_GN.update_func
%   Perform initial estimate of Homogeneous background
%    inv_model.inv_solve_abs_GN.do_starting_estimate  (default 1)
%
%   The signature for the update_func:
%    [img opt] = my_function(org, next, dx, fmin,res, opt)
%   where:
%    org  - the current estimate
%    next - the proposed next estimate as returned by line_optimize
%    dx   - the current derivative
%    fmin - the step size chosen by line_optimize
%    res  - the residual on "next" estimate
%    opt  - the option structure as described above
%    img  - the estimate to be used at next iteration
%
% See also: LINE_OPTIMIZE

% (C) 2010-2013 Andy Adler & Bartłomiej Grychtol.
% License: GPL version 2 or version 3
% $Id$


opt = parse_options(inv_model);
if opt.do_starting_estimate
    img = initial_estimate( inv_model, data1 );
else
    img = calc_jacobian_bkgnd( inv_model );
end

hp  = calc_hyperparameter( inv_model );
RtR = calc_RtR_prior( inv_model );
W   = calc_meas_icov( inv_model );
hp2RtR= hp*RtR;

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
  
  opt.line_optimize.hp2RtR = hp2RtR;
  [next fmin res] = line_optimize(img, dx, data1, opt.line_optimize);
  
  [img opt] = update_step(img, next, dx, fmin, res, opt);
  
  inv_model.jacobian_backgnd = img;
  RtR = calc_RtR_prior( inv_model );
  hp2RtR = hp*RtR;
end


function val = GN_objective_function(data0, data, img0, img, opt)
dv = calc_difference_data(data, data0, img0.fwd_model);
if ~isfield(img0, 'elem_data'), img0 = physics_data_mapper(img0); end
if ~isfield(img, 'elem_data'), img = physics_data_mapper(img); end
de = img0.elem_data - img.elem_data;
val = norm(dv) + 0.5 * de' * opt.hp2RtR * de;


function img = initial_estimate( imdl, data )
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
 
function [img opt] = update_step(org, next, dx, fmin,res, opt)
   if isfield(opt, 'update_func')
      [img opt] = feval(opt.update_func,org,next,dx,fmin,res,opt);
   else
      img = next;
   end
   
   
function opt = parse_options(imdl)

   opt.max_iter = 1;
   try
      opt.max_iter = imdl.parameters.max_iterations;
   end
   
   if isfield(imdl, 'inv_solve_abs_GN');
      opt = imdl.inv_solve_abs_GN;
   end
   
   if ~isfield(opt,'line_optimize')
      opt.line_optimize = [];
   end

   if isfield(opt, 'min_value')
      opt.line_optimize.min_value = opt.min_value;
   end
   
   if isfield(opt, 'max_value')
      opt.line_optimize.max_value = opt.max_value;
   end
   
   if ~isfield(opt, 'line_optimize') || ...
      ~isfield(opt.line_optimize, 'objective_func')
    % not sure this should be allowed to change
      opt.line_optimize.objective_func = ...
          @(d0,d,i0,i,o)GN_objective_function(d0,d,i0,i,o);
   end


