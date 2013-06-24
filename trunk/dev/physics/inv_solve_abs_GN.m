function img= inv_solve_abs_GN( inv_model, data0);
%INV_SOLVE_ABS_GN Absolute solver using Gauss Newton approximation
% img= inv_solve_abs_GN( inv_model, data0)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data0      => EIT data
%
% Parameters:
%   inv_model.parameters.max_iterations = N_max iter            (default 1)
%   inv_model.parameters.show_iterations  (print status lines)  (default 0)
% Parameters: (will override parameters field)
%   Maximum Iterations:
%    inv_model.inv_solve_abs_GN.max_iterations
%   Line Optimize function (more details below):
%    inv_model.inv_solve_abs_GN.line_optimize_func (default @line_optimize)
%   Line Optimize parameters:
%    inv_model.inv_solve_abs_GN.line_optimize.perturb            (optional)
%   Limits on solution:
%    inv_model.inv_solve_abs_GN.min_value                    (default -Inf)
%    inv_model.inv_solve_abs_GN.max_value                    (default +Inf)
%   Update function after each iteration (more details below):
%    inv_model.inv_solve_abs_GN.update_func                      (optional)
%   Perform initial estimate (more details below):
%    inv_model.inv_solve_abs_GN.do_starting_estimate            (default 1)
%
%   Signature for line_optimize_func
%    [next fmin res] = my_line_optimize(org, dx, data0, opt)
%   where:
%    next - the image to use in the next iteration
%    fmin - step size along dx that minimizes the objective function
%    res  - the value of the objective function for the next image
%    org  - the current estimate
%    dx   - the current derivative (step direction)
%    opt  - the options structure specified in 
%           inv_model.inv_solve_abs_GN.line_optimize (can be used to pass
%           additonal arguments, stores the GN objective function in 
%           opt.objective_func)
%   Note that changing the line optimize function can substantially alter
%   the behavior of the solver, especially if a different objective
%   function is used.
%
%   Signature for update_func:
%    [img opt] = my_update_func(org, next, dx, fmin,res, opt)
%   where:
%    org  - the current estimate
%    next - the proposed next estimate as returned by line_optimize
%    dx   - the current derivative (step direction)
%    fmin - the step size chosen by line_optimize
%    res  - the residual on "next" estimate
%    opt  - the option structure as described above
%    img  - the estimate to be used at next iteration
%
% By default, the starting estimate will be based on
% inv_model.jacobian_bkgnd, but scaled such as to best fit data0. 
% Set opt.do_initial_estimate to 0 to use inv_model.jacobian_bkgnd without
% modification. 
% The initial estimate serves as a prior for the reconstruction, deviation
% from which is penalized.
%
% See also: LINE_OPTIMIZE

% (C) 2010-2013 Andy Adler & Bartłomiej Grychtol.
% License: GPL version 2 or version 3
% $Id$


opt = parse_options(inv_model);
if opt.do_starting_estimate
    img = initial_estimate( inv_model, data0 );
else
    img = calc_jacobian_bkgnd( inv_model );
end

img = physics_param_mapper(img);

hp  = calc_hyperparameter( inv_model );
RtR = calc_RtR_prior( inv_model );
W   = calc_meas_icov( inv_model );
hp2RtR= hp*RtR;

img0 = physics_param_mapper(img);
opt.line_optimize.meas_icov = calc_meas_icov( inv_model);
for i = 1:opt.max_iter
  vsim = fwd_solve( img ); 
  dv = calc_difference_data( vsim , data0, img.fwd_model);
  J = calc_jacobian( img );

  % create param
  tmp = physics_param_mapper(img);
  RDx = hp2RtR*(img0.params - tmp.params);
  dx = (J'*W*J + hp2RtR)\(J'*dv + RDx);
  
  opt.line_optimize.hp2RtR = hp2RtR;
  [next, fmin, res] = ...
      feval(opt.line_optimize_func,img, dx, data0, opt.line_optimize);

  if opt.show_iterations
     eidors_msg('#%02d residual=%.3g', i, res, 1);
  end
  
  [img, opt] = update_step(img, next, dx, fmin, res, opt);
  
  inv_model.jacobian_backgnd = physics_param_mapper(img,1);
  RtR = calc_RtR_prior( inv_model );
  hp2RtR = hp*RtR;
end


function val = GN_objective_function(data0, data, img0, img, opt)
   dv = calc_difference_data(data, data0, img0.fwd_model);
   if ~isfield(img0, 'params'), img0 = physics_param_mapper(img0); end
   if ~isfield(img, 'params'), img = physics_param_mapper(img); end
   de = img0.params - img.params;
   val = 0.5*( dv'*opt.meas_icov*dv + de' * opt.hp2RtR * de);


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
   img = physics_param_mapper(img);
   
   if isfield(img.fwd_model,'coarse2fine');
      % TODO: the whole coarse2fine needs work here.
      %   what happens if c2f doesn't cover the whole region
      
      % TODO: this assumes that params is elem_data...
      
      % TODO: the two cases are very different. c2f case should match other
      nc = size(img.fwd_model.coarse2fine,2);
      img.params = mean(img.params)*ones(nc,1)*pf(1);
      
      warning('c2f needs work');
   else
      img.params = img.params*pf(1);
   end
   
   % remove elem_data
   img = physics_param_mapper(img,1);
 
function [img opt] = update_step(org, next, dx, fmin,res, opt)
   if isfield(opt, 'update_func')
      [img opt] = feval(opt.update_func,org,next,dx,fmin,res,opt);
   else
      img = next;
   end
   
   
function opt = parse_options(imdl)
   try
       % for any general options
       opt = imdl.parameters;
   end
   
   opt.max_iter = 1;

   try
      opt.max_iter = imdl.parameters.max_iterations;
   end
   
   if isfield(imdl, 'inv_solve_abs_GN');
      fnames = fieldnames(imdl.inv_solve_abs_GN);
      for i = 1:length(fnames)
          opt.(fnames{i}) = imdl.inv_solve_abs_GN.(fnames{i});
      end
   end
   
   if ~isfield(opt,'line_optimize_func')
      opt.line_optimize_func = @line_optimize;
   end
   
   if ~isfield(opt,'line_optimize')
      opt.line_optimize = [];
   end

   if ~isfield(opt,'show_iterations')
      opt.show_iterations = 0;
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
      opt.line_optimize.objective_func = @GN_objective_function;       
   end
   
   if ~isfield(opt,'do_starting_estimate')
       opt.do_starting_estimate = 1;
   end

