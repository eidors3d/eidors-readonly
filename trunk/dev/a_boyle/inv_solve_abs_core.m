function img= inv_solve_abs_core( inv_model, data0);
%INV_SOLVE_ABS_CORE Absolute solver using a generic iterative algorithm
%
% NOTE that, in general, this function should not be called directly.
%
% See also: INV_SOLVE_ABS_GN, INV_SOLVE_ABS_GN_LOGC,
%           INV_SOLVE_ABS_CG, INV_SOLVE_ABS_CG_LOGC,
%           LINE_SEARCH_O2, LINE_SEARCH_ONM2
%
% (C) 2010-2013 Andy Adler & Bartłomiej Grychtol, Nolwenn Lespare, Alistair Boyle.
% License: GPL version 2 or version 3
% $Id$

% TODO this documentation is way, way out of date... it needs a serious rewrite!
%
% img= inv_solve_abs_GN( inv_model, data0)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data0      => EIT data
%
% Parameters:
%   inv_model.parameters.max_iterations = N_max iter            (default 1)
%   inv_model.parameters.show_iterations  (print status lines)  (default 0)
%   Line Search function (more details below):
%    inv_model.parameters.line_search_func        (default @line_search_o2)
%   Line Search parameters:
%    inv_model.parameters.line_search_perturb                    (optional)
%   Limits on solution:
%    inv_model.parameters.min_value                          (default -Inf)
%    inv_model.parameters.max_value                          (default +Inf)
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

%--------------------------
% UNIT_TEST?
if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); img = do_unit_test; return; end

%--------------------------
opt = parse_options(inv_model);
%if opt.do_starting_estimate
%    img = initial_estimate( inv_model, data0 );
%else
[inv_model, opt] = append_c2f_background(inv_model, opt);
img = calc_jacobian_bkgnd( inv_model );
% TODO does calc_jacobian_bkgnd ignore 'physics' right now.. that might screw things up pretty good!
img = physics_data_mapper(img); % move data from whatever 'physics' to img.elem_data


% mk_image doesn't handle the c2f
% TODO move this from here to mk_image so we get an img with the correct number of elem_data
if isfield(inv_model, 'fwd_model') && ...
   isfield(inv_model.fwd_model, 'coarse2fine') && ...
   length(img.elem_data) ~= size(inv_model.fwd_model.coarse2fine,2)

   c2f = inv_model.fwd_model.coarse2fine;
   img.elem_data = c2f \ img.elem_data;
   bg = mean(img.elem_data);
   %img.elem_data = ones(size(c2f,2),1)*bg;
   if opt.verbose > 1
      fprintf('  c2f: correcting mk_image elem_data size %d -> %d (av %0.1f Ohm.m)\n', size(c2f), 1/bg);
      disp('    TODO this fix should be moved to mk_image()');
   end
end


% precalculate some of our matrices if required
hp = init_hp(inv_model, opt);
W  = init_meas_icov(inv_model, opt);
N = init_normalization(inv_model.fwd_model, opt);

% map data and measurements to working types
%  convert elem_data
img.elem_data = map_data(img.elem_data, img.current_physics, opt.elem_working);
img.current_physics = opt.elem_working;
%  convert measurement data
if ~isstruct(data0)
  d = data0;
  data0 = struct;
  data0.meas = d;
  data0.type = 'data';
end
data0.meas = map_meas(data0.meas, N, opt.meas_input, opt.meas_working);
data0.current_physics = opt.meas_working;

% now get on with 
img0 = img;
RtR = 0; k = 0; dv = []; de = []; sx = 0; r = 0; stop = 0; % general init
residuals = zeros(opt.max_iterations,3); fig_r = []; % for residuals plots
while ~stop
  if opt.verbose > 1
     fprintf('  iteration %d\n', k+1)
  end

  % calculate the Jacobian, possibly update RtR as well if required
  J = update_jacobian(img, opt);
  RtR = update_RtR(RtR, inv_model, k, img, opt);

  % update change in element data from the prior de and
  % the measurement error dv
  [dv, opt] = update_dv(dv, img, img0, data0, N, opt);
  de = update_de(de, img, img0, opt);

  % now find the residual, quit if we're done
  [stop, k, r, fig_r] = update_residual(dv, de, W, hp, RtR, k, r, fig_r, opt);
  if stop
     break;
  end

  % determine the next search direction sx
  %  dx is specific to the algorithm, generally "downhill"
  dx = update_dx(J, W, hp, RtR, dv, de, opt);
  % choose beta, beta=0 unless doing Conjugate Gradient
  beta = update_beta(dx, opt);
  % sx_k = dx_k + beta * sx_{k-1}
  sx = update_sx(dx, beta, sx);

  % line search for alpha, leaving the final selection as img
  [alpha, img, dv, opt] = update_alpha(img, sx, data0, img0, N, W, hp, RtR, dv, opt);
  % fix max/min values for x, clears dx if limits are hit, where
  % a cleared dv will trigger a recalculation of dv at the next update_dv()
  [img, dv] = update_img_using_limits(img, img0, data0, N, dv, opt);

  if opt.verbose > 5
     show_fem_iter(k, img, inv_model, opt);
  end
end
[img, opt] = strip_c2f_background(img, opt);
% check we're returning the right size of data
if isfield(inv_model, 'rec_model')
  img.fwd_model = inv_model.rec_model;
end
if opt.verbose > 1
   fprintf('  %d fwd_solves required for this solution in %d iterations\n', ...
           opt.fwd_solutions, k);
end
if size(img.elem_data,1) ~= size(img.fwd_model.elems,1)
   error('img data (%d) does not match FEM mesh elements (%d)', ...
         length(img.elem_data), length(img.fwd_model.elems));
end
% convert data for output
img.elem_data = map_data(img.elem_data, img.current_physics, opt.elem_output);
img.current_physics = opt.elem_output;
%img.physics_data_mapper = opt.elem_output; %TODO which is it: current_physics or physics_data_mapper???
%img = physics_data_mapper(img, 1); % move data from img.elem_data to whatever 'physics'

function W = init_meas_icov(inv_model, opt)
   W = 1;
   if opt.calc_meas_icov
      if opt.verbose > 1
         disp('  calc measurement inverse covariance W');
      end
      W   = calc_meas_icov( inv_model );
   end

function hp = init_hp(inv_model, opt)
   hp = 0;
   if opt.calc_hyperparameter
      if opt.verbose > 1
         disp('  calc regularization hyperparameter(s)');
      end
      hp  = calc_hyperparameter( inv_model );
   end

function N = init_normalization(fmdl, opt)
   % precalculate the normalization of the data if required (apparent resistivity)
   N = 1;
   if opt.normalize_data
      if opt.verbose > 1
         disp('  calc measurement normalization matrix N');
      end
      N = feval(opt.normalize_data_func, fmdl);
   end

% r_km1: previous residual, if its the first iteration r_km1 = inf
% r_k: new residual
% fig_r: the handle to the residual plot if used
function [stop, k, r, fig_r] = update_residual(dv, de, W, hp, RtR, k, r, fig_r, opt)
  stop = 0;
  % update iteration count
  k = k+1;

  % update residual estimate
  if k == 1
     r = zeros(opt.max_iterations, 3);
     r_km1 = inf;
  else
     r_km1 = r(k-1, 1);
  end
  r_k = feval(opt.residual_func, dv, de, W, hp, RtR);
  % save residual for next iteration
  r(k,1) = r_k;

  % now do something with that information
  if opt.verbose > 1
    fprintf('    calc residual\n');
    fprintf('      r[%d] = %0.3e\n', k, r_k);
    dr = (r_k - r_km1)/r_k;
    fprintf('      dr (r[%d]-r[%d])/r[%d] = %0.2f (%0.1e)\n', k, k-1, k, dr, dr);
  end
  if opt.plot_residuals
     %         optimization_criteria, data misfit, roughness
     r(k,2:3) = [(dv'*dv)/2 (de'*de)/2];
     if isempty(fig_r)
        fig_r = figure();
     else
        figure(fig_r);
     end
     x = 1:k;
     y = r(x, :);
     y = y ./ repmat(max(y,[],1),size(y,1),1) * 100;
     plot(x, y, 'o-', 'linewidth', 2, 'MarkerSize', 10);
     title('residuals');
     axis tight;
     ylabel('residual (% of max)');
     xlabel('iteration');
     legend('residual','meas. misfit','prior misfit');
     legend('Location', 'EastOutside');
     drawnow;
  end
  
  % evaluate stopping criteria
  if k > opt.max_iterations
     stop = 1;
  end
  % TODO return 'measurement residual' & 'roughness' for progress plot, as well
  if r_k < opt.tol + opt.ntol
     if opt.verbose > 1
        fprintf('  terminated at iteration %d\n',k);
        fprintf('    residual tolerance (%0.1e) achieved\n', opt.tol + opt.ntol);
     end
     stop = 1;
  end
  if (k > opt.dtol_iter) && ((r_k - r_km1)/r_k > opt.dtol + 2*opt.ntol)
     if opt.verbose > 1
        fprintf('  terminated at iteration %d (iterations not improving)\n', k);
        fprintf('    residual slope tolerance (%0.1e) exceeded\n', opt.dtol + 2*opt.ntol);
     end
     stop = 1;
  end

% for Conjugate Gradient, else beta = 0
function beta = update_beta(dx, opt);
   beta = 0;

% update the search direction
% for Gauss-Newton
%   sx_k = dx_k
% for Conjugate-Gradient
%   sx_k = dx_k + beta * sx_{k-1}
function sx = update_sx(dx, beta, sx);
   sx = dx + beta * sx;

function RtR = update_RtR(RtR, inv_model, k, img, opt)
   % TODO sometimes (with Noser?) this requires the Jacobian, could this be done more efficiently?
   % add a test function to determine if img.elem_data affects RtR, skip this if independant
   % TODO we could detect in the opt_parsing whether the calc_RtR_prior depends on 'x' and skip this if no
   if opt.calc_RtR_prior
      if opt.verbose > 1
         try RtR_str = func2str(inv_model.RtR_prior);
         catch RtR_str = inv_model.RtR_prior;
         end
         fprintf('    calc regularization RtR (%s)\n', RtR_str);
      end
      RtR = calc_RtR_prior_wrapper(inv_model, img, opt);
   else
      RtR = 0;
   end

function J = update_jacobian(img, opt)
   if opt.calc_jacobian % TODO skipping the Jacobian isn't really an option ever... delete this 'if'
      img.elem_data = map_data(img.elem_data, img.current_physics, 'conductivity');
      img.current_physics = 'conductivity';
      if(opt.verbose > 1)
         try J_str = func2str(img.fwd_model.jacobian);
         catch J_str = img.fwd_model.jacobian;
         end
         fprintf('    calc Jacobian J @ current x (%s)\n', J_str);
      end
      % scaling if we are working in something other than direct conductivity
      S = feval(opt.calc_jacobian_scaling_func, img.elem_data); % chain rule
      % finalize the jacobian
      J = calc_jacobian( img ) * S;
   else
      J = [];
   end
% -------------------------------------------------
% Chain Rule Products for Jacobian Translations
% x is conductivity, we want the chain rule to translate the
% Jacobian of conductivity to conductivity on resistivity or
% logs of either.
% This chain rule works out to a constant.
%
% d log_b(x)     1          d x
% ---------- = ------- , ---------- = x ln(b)
%     d x      x ln(b)   d log_b(x)
function S = dx_dlogx(x);
   S = diag(x);
function S = dx_dlog10x(x);
   S = diag(x) * log(10);
% resistivity
% dx      d x   -1
% ----- = --- = ---, y = 1/x --> -(x^2)
% d 1/x   d y   y^2
function S = dx_dy(x);
   S = diag(-(x.^2));
% then build the log versions of conductivity by combining chain rule products
function S = dx_dlogy(x);
%   S = dx_dy(x) * dy_dlogy(x);
%     = -(x^2) * 1/x = -x
   S = diag(-x);
function S = dx_dlog10y(x);
%   S = dx_dy(x) * dy_dlog10y(x);
%     = -(x^2) * 1/(ln(10) x) = -x / ln(10)
   S = diag(-x/log(10));
% ... some renaming to make things understandable above: x = 1/y
%function S = dy_dlogy(x);
%   S = dx_dlogx(1./x);
%function S = dy_dlog10y(x);
%   S = dx_dlog10x(1./x);
% -------------------------------------------------


function [alpha, img, dv, opt] = update_alpha(img, sx, data0, img0, N, W, hp, RtR, dv, opt)
  if(opt.verbose > 1)
     disp('    line search');
  end
  [alpha, img, dv, opt] = feval(opt.line_search_func, img, sx, data0, img0, N, W, hp, RtR, dv, opt);

function [img, dv] = update_img_using_limits(img, img0, data0, N, dv, opt)
  % fix max/min values for x
  if opt.max_value ~= +inf
     lih = find(img.elem_data > opt.max_value);
     img.elem_data(lih) = opt.max_value;
     if opt.verbose > 1
        fprintf('    limit max(x)=%e for %d elements\n', opt.max_value, length(lih));  
     end
     dv = []; % dv is now invalid since we changed the conductivity
  end
  if opt.min_value ~= -inf
     lil = find(img.elem_data < opt.min_value);
     img.elem_data(lil) = opt.min_value;
     if opt.verbose > 1
        fprintf('    limit min(x)=%e for %d elements\n', opt.min_value, length(lil));  
     end
     dv = []; % dv is now invalid since we changed the conductivity
  end
  % update voltage change estimate if the limit operation changed the img data
  [dv, opt] = update_dv(dv, img, img0, data0, N, opt, '(dv out-of-date)');

function  de = update_de(de, img, img0, opt)
   % probably not the most robust check for whether this is the first update
   % but this ensures that we get exactly zero for the first iteration and not
   % a set of values that has numeric floating point errors that are nearly zero
   if isempty(de) % first iteration
      % data hasn't changed yet!
      de = zeros(size(img0.elem_data));
   else
      de = img0.elem_data - img.elem_data;
   end

function [dv, opt] = update_dv(dv, img, img0, data0, N, opt, reason)
   % estimate current error as a residual
   if ~isempty(dv) % need to calculate dv...
      return;
   end
   if nargin < 7
      reason = '';
   end
   if opt.verbose > 1
      disp(['    fwd_solve ', reason]);
   end
   [dv, opt] = update_dv_core(img, img0, data0, N, opt);

% also used by the line search as opt.line_search_dv_func
function [dv, opt] = update_dv_core(img, img0, data0, N, opt)
   img.elem_data = map_data(img.elem_data, img.current_physics, 'conductivity');
   img.current_physics = 'conductivity';
   data = fwd_solve(img);
   try   current_meas_physics = data.current_physics;
   catch current_meas_physics = 'voltage';
   end
   data.meas = map_meas(data.meas, N, current_meas_physics, opt.meas_working);
   data.current_physics = opt.meas_working;
   opt.fwd_solutions = opt.fwd_solutions +1;
   dv = calc_difference_data(data, data0, img0.fwd_model);

function show_fem_iter(k, img, inv_model, opt)
  [img, opt] = strip_c2f_background(img, opt);
  % check we're returning the right size of data
  if isfield(inv_model, 'rec_model')
    img.fwd_model = inv_model.rec_model;
  end
  bg = 1;
  img.calc_colours.ref_level = bg;
  img.calc_colours.clim = bg;
  img.calc_colours.cb_shrink_move = [0.3,0.6,0.02]; % move color bars
  figure; show_fem(img, 1);
  title(sprintf('iter=%d',k));

% TODO confirm that GN line_search_onm2 is using this residual calculation (preferably, directly)
function residual = GN_residual(dv, de, W, hp, RtR)
   % we operate on whatever the iterations operate on (log data, resistance, etc) + perturb(i)*dx
   hp2RtR = hp*RtR;
   residual = 0.5*( dv'*W*dv + de'*hp2RtR*de);

function RtR = calc_RtR_prior_wrapper(inv_model, img, opt)
   inv_model.jacobian_backgnd = img;
   RtR = calc_RtR_prior( inv_model );
   if size(RtR,1) < length(img.elem_data)
     ne = length(img.elem_data) - size(RtR,1);
     RtR(end+1:end+ne, end+1:end+ne) = 1;
     if opt.verbose > 1
        fprintf('    c2f: adjusting RtR by appending %d rows/cols\n', ne);
        disp('      TODO move this fix, or something like it to calc_RtR_prior -- this fix is a quick HACK to get things to run...');
     end
   end

%function img = initial_estimate( imdl, data )
%   img = calc_jacobian_bkgnd( imdl );
%   vs = fwd_solve(img);
%
%   if isstruct(data)
%      data = data.meas;
%   else
%     meas_select = [];
%     try
%        meas_select = imdl.fwd_model.meas_select;
%     end
%     if length(data) == length(meas_select)
%        data = data(meas_select);
%     end
%   end
%
%   pf = polyfit(data,vs.meas,1);
%
%   % create elem_data
%   img = physics_data_mapper(img);
%
%   if isfield(img.fwd_model,'coarse2fine');
%      % TODO: the whole coarse2fine needs work here.
%      %   what happens if c2f doesn't cover the whole region
%
%      % TODO: the two cases are very different. c2f case should match other
%      nc = size(img.fwd_model.coarse2fine,2);
%      img.elem_data = mean(img.elem_data)*ones(nc,1)*pf(1);
%   else
%      img.elem_data = img.elem_data*pf(1);
%   end
%
%   % remove elem_data
%%   img = physics_data_mapper(img,1);
%
%function [img opt] = update_step(org, next, dx, fmin,res, opt)
%   if isfield(opt, 'update_func')
%      [img opt] = feval(opt.update_func,org,next,dx,fmin,res,opt);
%   else
%      img = next;
%   end


function opt = parse_options(imdl)
   try
      % for any general options
      opt = imdl.parameters;
   catch
      opt = struct;
   end

   % verbosity, debug output
   % 0: quiet
   % 1: print iteration count
   % 2: print details as the algorithm progresses
   if ~isfield(opt,'verbose')
      opt.verbose = 4;
   end
   if(opt.verbose > 1)
      fprintf('  setting default parameters\n');
   end
   % we track how many fwd_solves we do since they are the most expensive part of the iterations
   opt.fwd_solutions = 0;

   if ~isfield(opt, 'residual_func') % the objective function
      opt.residual_func = @GN_residual; % r = f(dv, de, W, hp, RtR)
   end

   % calculation of update components
   if ~isfield(opt, 'update_func')
      opt.update_func = @GN_update; % dx = f(J, W, hp, RtR, dv, de)
   end
   % figure out if things need to be calculated
   if ~isfield(opt, 'calc_jacobian') % derivative of the objective function
      opt.calc_jacobian = 0; % J
   end
   if ~isfield(opt, 'calc_meas_icov') % derivative of the objective function
      opt.calc_meas_icov = 0; % W
   end
   if ~isfield(opt, 'calc_RtR_prior') % derivative of the objective function
      opt.calc_RtR_prior = 0; % RtR
   end
   if ~isfield(opt, 'calc_hyperparameter') % derivative of the objective function
      opt.calc_hyperparameter = 0; % hp
   end
%   try
      if opt.verbose > 1
         fprintf('    examining function %s(...) for required arguments\n', func2str(opt.update_func));
      end
      % ensure that necessary components are calculated
      % opt.update_func: dx = f(J, W, hp, RtR, dv, de)
      args = function_depends_upon(opt.update_func, 6);
      if args(1) == 1
         opt.calc_jacobian = 1;
      end
      if args(2) == 1
         opt.calc_meas_icov = 1;
      end
      if args(3) == 1
         opt.calc_hyperparameter = 1;
      end
      if args(4) == 1
         opt.calc_RtR_prior = 1;
      end
%   catch
%      error('exploration of function %s via function_depends_upon() failed', func2str(opt.update_func));
%   end
      

   % stopping criteria, solution limits
   if ~isfield(opt, 'max_iterations')
      opt.max_iterations = 10;
   end
   if ~isfield(opt, 'ntol')
      opt.ntol = eps; % attempt to quantify numeric machine precision
   end
   if ~isfield(opt, 'tol')
      opt.tol = 0; % terminate iterations if residual is less than tol
   end
   if ~isfield(opt, 'dtol')
      % terminate iterations if residual slope is greater than dtol
      % generally, we would want dtol to be -0.01 (1% decrease) or something similar
      %  ... as progress levels out, stop working
      %opt.dtol = +inf;
      opt.dtol = -1e-4; % --> -0.01% slope (really slow)
   end
   if ~isfield(opt, 'dtol_iter')
      %opt.dtol_iter = inf; % ignore dtol for dtol_iter iterations
      opt.dtol_iter = 0; % use dtol from the begining
   end
   if ~isfield(opt, 'min_value')
      opt.min_value = -inf; % min elem_data value
   end
   if ~isfield(opt, 'max_value')
      opt.max_value = +inf; % max elem_data value
   end
   % provide a graphical display of the line search values & fit
   if ~isfield(opt, 'plot_residuals')
      if opt.verbose > 2
         opt.plot_residuals = 1;
      else
         opt.plot_residuals = 0;
      end
   end
   if opt.plot_residuals ~= 0
      disp('  residual plot (updated per iteration) are enabled, to disable them set');
      disp('    inv_model.parameters.plot_residuals=0');
   end

   % line search
   if ~isfield(opt,'line_search_func')
      %opt.line_search_func = @line_search_o2; % img = f(img, dx, opt) TODO out of date!!!
      opt.line_search_func = @line_search_onm2; % img = f(img, dx, opt) TODO out of date!!!
   end
   if ~isfield(opt,'line_search_dv_func')
      opt.line_search_dv_func = @update_dv_core;
      % [dv, opt] = update_dv_core(img, img0, data0, N, opt)
   end
   % an initial guess for the line search step sizes, may be modified by line search
   if ~isfield(opt,'line_search_args') || ...
      ~isfield(opt.line_search_args, 'perturb')
      fmin = 1/4; % arbitrary starting guess
      %opt.line_search_args.perturb = [0 fmin/4 fmin/2 fmin fmin*2 fmin*4];
      %opt.line_search_args.perturb = [0 fmin/4 fmin fmin*4];
      %opt.line_search_args.perturb = [0 0.1 0.5 0.7 1.0];
      opt.line_search_args.perturb = [0 0.1 0.7 0.9 1.0];
      %opt.line_search_args.perturb = [0 0.1 0.9 1.0];
   end
   % provide a graphical display of the line search values & fit
   if ~isfield(opt,'line_search_args') || ...
      ~isfield(opt.line_search_args, 'plot')
      if opt.verbose > 3
         opt.line_search_args.plot = 1;
      else
         opt.line_search_args.plot = 0;
      end
   end
   if opt.line_search_args.plot ~= 0
      disp('  line search plots (per iteration) are enabled, to disable them set');
      disp('    inv_model.parameters.line_search_args.plot=0');
   end

   % fixed elem_data, ignored by the algorithm
   if ~isfield(opt, 'ignore_elems')
      opt.ignore_elems = [];
   end

   % background
   % if > 0, this is the elem_data that holds the background
   % this is stripped off just before the iterations complete
   if ~isfield(opt, 'c2f_background')
     if isfield(imdl, 'fwd_model') && isfield(imdl.fwd_model, 'coarse2fine')
        opt.c2f_background = -1; % possible: check if its required later
     else
        opt.c2f_background = 0;
     end
   end

   % DATA CONVERSION settings
   % elem type for the initial estimate is based on calc_jacobian_bkgnd which returns an img
   if ~isfield(opt, 'elem_working')
      opt.elem_working = 'conductivity';
   end
   if ~isfield(opt, 'elem_output')
      opt.elem_output = 'conductivity';
   end
   if ~isfield(opt, 'meas_input')
      opt.meas_input = 'voltage';
   end
   if ~isfield(opt, 'meas_working')
      opt.meas_working = 'voltage';
   end

   % JACOBIAN CHAIN RULE conductivity -> whatever
   % where x = conductivity at this iteration
   %       S = a scaling matrix, generally a diagonal matrix of size matching Jacobian columns
   % Jn = J * S;
   % if not provided, determine based on 'elem_working' type
   if ~isfield(opt, 'calc_jacobian_scaling_func')
      switch opt.elem_working
         case 'conductivity'
            opt.calc_jacobian_scaling_func = @ret1_func;  % S = f(x)
         case 'log_conductivity'
            opt.calc_jacobian_scaling_func = @dx_dlogx;   % S = f(x)
         case 'log10_conductivity'
            opt.calc_jacobian_scaling_func = @dx_dlog10x; % S = f(x)
         case 'resistivity'
            opt.calc_jacobian_scaling_func = @dx_dy;      % S = f(x)
         case 'log_resistivity'
            opt.calc_jacobian_scaling_func = @dx_dlogy;   % S = f(x)
         case 'log10_resistivity'
            opt.calc_jacobian_scaling_func = @dx_dlog10y; % S = f(x)
         otherwise
            warning('unrecognized opt.elem_working type, please provide a opt.calc_jacobian_scaling_func -- assuming no scalling and continuing');
            opt.calc_jacobian_scaling_func = @ret1_func; % S = f(x)
      end
   end

   % input handling -> conversion of measurements (data0)
   if ~isfield(opt, 'normalize_data_func') % how do we normalize data? use this function
      if strcmp(opt.meas_input,   'apparent_resistivity') || ...
         strcmp(opt.meas_working, 'apparent_resistivity')
         opt.normalize_data_func = @calc_normalization_apparent_resistivity; % N = f(fmdl)
      else
         opt.normalize_data_func = @ret1_func; % N = f(fmdl)
      end
   end
   if isfield(opt, 'normalize_data_func')
      opt.normalize_data = 1;
   else
      opt.normalize_data = 0;
   end

function check_matrix_sizes(J, W, hp, RtR, dv, de, opt)
   % assuming our equation looks something like
   % dx = (J'*W*J + hp2RtR)\(J'*dv + hp2RtR*de);
   % check that all the matrix sizes are correct
   ne = size(de,1);
   nv = size(dv,1);
   if size(de,2) ~= 1
      error('de cols (%d) not equal 1', size(de,2));
   end
   if size(dv,2) ~= 1
      error('dv cols (%d) not equal 1', size(dv,2));
   end
   if opt.calc_meas_icov && ...
      any(size(W) ~= [nv nv])
      error('W size (%d rows, %d cols) is incorrect (%d rows, %d cols)', size(W), nv, nv);
   end
   if opt.calc_jacobian && ...
      any(size(J) ~= [nv ne])
      error('J size (%d rows, %d cols) is incorrect (%d rows, %d cols)', size(J), nv, ne);
   end
   if opt.calc_RtR_prior && ...
      any(size(RtR) ~= [ne ne])
      error('RtR size (%d rows, %d cols) is incorrect (%d rows, %d cols)', size(RtR), ne, ne);
   end

function N = calc_normalization_apparent_resistivity(fmdl)
  img1 = mk_image(fmdl,1);
  vh1  = fwd_solve(img1);
  N    = diag(1./vh1.meas);

%  normalisation= 1./vh1.meas;
%  N= speye(length(normalisation));
%  N(1:size(N,1)+1:size(N,1)*size(N,1))= normalisation;
%  sz_N = size(N)
%  figure; spy(N)

function dx = update_dx(J, W, hp, RtR, dv, de, opt)
   if(opt.verbose > 1)
      fprintf( '    calc step size dx');
   end

   % TODO move this outside the inner loop of the iterations, it only needs to be done once
   check_matrix_sizes(J, W, hp, RtR, dv, de, opt)

   de(opt.ignore_elems) = 0;
   if ~isempty(RtR)
      RtR(:,opt.ignore_elems) = 0;
   end

   % do the update step direction calculation
   dx = feval(opt.update_func, J, W, hp, RtR, dv, de);
   % ignore any fixed value elements
   dx(opt.ignore_elems) = 0;

   if(opt.verbose > 1)
      fprintf(', ||dx||=%0.2e\n', norm(dx));
   end

function dx = GN_update(J, W, hp, RtR, dv, de)
   hp2RtR = hp*RtR;
   % the actual update
   dx = (J'*W*J + hp2RtR)\(J'*dv + hp2RtR*de);

% for each argument, returns 1 if the function depends on it, 0 otherwise
% 'zero' arguments do not need to be calculated since they don't get used
function args = function_depends_upon(func, argn)
   % build function call
   str = sprintf('%s(',func2str(func));
   args = zeros(argn,1);
   for i = 1:argn-1
      str = [str sprintf('a(%d),',i)];
   end
   str = [str sprintf('a(%d))',argn)];
   % baseline
   a = ones(argn,1)*2;
   x = eval(str);
   % now check for a difference at each argument
   for i = 1:argn
      a = ones(argn,1)*2;
      a(i) = 0;
      y = eval(str);
      if any(x ~= y)
         args(i) = 1;
      end
   end

% this function just passes data from its input to its output
function out = null_func(in, arg1);
   out = in;

% this function always returns one
function out = ret1_func(arg1, arg2);
   out = 1;

% if required, expand the coarse-to-fine matrix to cover the background of the image
% this is removed at the end of the iterations
function [inv_model, opt] = append_c2f_background(inv_model, opt)
    % either there is already a background
    % or none is required --> -1 means we go to work building one
    if opt.c2f_background >= 0
      return
    end
    % check that all elements get assigned a conductivity
    % through the c2f conversion
    c2f = inv_model.fwd_model.coarse2fine; % coarse-to-fine mesh mapping
    nf = size(inv_model.fwd_model.elems,1); % number of fine elements
    nc = size(c2f,2); % number of coarse elements
    % ... each element of fel aught to sum to '1' since the
    % elem_data is being assigned from a continuous unit
    % surface value
    % now, find any fine elements that are not fully
    % mapped between the two meshes (<1) w/in a tolerance
    % related to the number of additions in the summation
    fel = sum(c2f,2); % collapse mapping onto the fwd_model elements
    n = find(fel < 1 - (1e3+nc)*eps);
    % ... 1e3 is a fudge factor since we don't care too much
    %     about small area mapping errors
    % if we do have some unassigned elements,
    % expand c2f and add a background element to the 'elem_data'
    if length(n) ~= 0
      if(opt.verbose > 1)
        fprintf('  c2f: adding background conductivity to %d\n    fwd_model elements not covered by rec_model\n', length(n));
      end
      c2f(n,nc+1) = 1 - fel(n);
      inv_model.fwd_model.coarse2fine = c2f;
      opt.c2f_background = nc+1;
    end

function [img, opt] = strip_c2f_background(img, opt)
    if opt.c2f_background <= 0
      return;
    end
    e = opt.c2f_background;
    img.elem_data_background = img.elem_data(e);
    img.elem_data(e) = [];
    img.fwd_model.coarse2fine(:,e) = [];
    opt.c2f_background = 0;
    if(opt.verbose > 1)
      fprintf('  background conductivity: %0.1f Ohm.m\n', 1./img.elem_data_background);
    end

function img = fix_c2f_calc_jacobian_backgnd(img, opt)
    % finally, we need to fix the img.elem_data returned by
    % mk_image/calc_jacobian_bkgnd, since it doesn't know
    % about coarse2fine.. yet TODO
    % somewhat arbitrarily, we choose the background value
    % as the mean background value.. we can do better TODO
    nc = size(c2f,2);
    if nc ~= length(img.elem_data)
      if(opt.verbose > 1)
        fprintf('  c2f: correcting img.elem_data to coarse size %d -> $d as mean value\n', length(img.elem_data), nc);
        fprintf('  c2f: TODO calc_jacobian_backgnd() should handle c2f!\n');
      end
      bg = mean(img.elem_data);
      img.elem_data = bg*ones(nc,1);
    end

function b = has_physics(s)
b = false;
if isstruct(s)
   b = any(ismember(fieldnames(s),supported_physics));
end

function x = map_data(x, in, out)
   if strcmp(in, out) % in == out
      return; % do nothing
   end

   % resistivity to conductivity conversion
   % we can't get here if in == out
   % we've already checked for log convserions on input or output
   if (strcmp(in,  'resistivity') || strcmp(in,  'conductivity')) && ...
      (strcmp(out, 'resistivity') || strcmp(out, 'conductivity'))
      x = 1./x; % conductivity <-> resistivity
   % log conversion
   elseif ~strcmp(in(1:6), out(1:6))
      % log_10 x -> x
      if strcmp(in(1:6), 'log10_')
         x = map_data(10.^x, in(7:end), out);
      % ln x -> x
      elseif strcmp(in(1:4), 'log_')
         x = map_data(exp(x), in(5:end), out);
      % x -> log_10 x
      elseif strcmp(out(1:6), 'log10_')
         x = log10(map_data(x, in, out(7:end)));
      % x -> ln x
      elseif strcmp(out(1:4), 'log_')
         x = log(map_data(x, in, out(5:end)));
      else
         error(sprintf('unknown conversion (log conversion?) %s - > %s', in, out));
      end
   else
      error('unknown conversion %s -> %s', in, out);
   end
   if any(isinf(x)) || any(isnan(x))
      warning('conversion results in Inf or NaN values');
   end

function b = map_meas(b, N, in, out)
   if strcmp(in, out) % in == out
      return; % do nothing
   end

   % log conversion
   if ~strcmp(in(1:4), out(1:4))
      % log -> natural
      if strcmp(in(1:4), 'log_')
         b = map_meas(10.^b, in(5:end), out);
         return;
      % natural -> log
      elseif strcmp(out(1:4), 'log_')
         b = log10(map_meas(b, in, out(5:end)));
         return;
      end
   end

   % resistivity to conductivity conversion
   % we can't get here if in == out
   if     strcmp(in, 'voltage') && strcmp(out, 'apparent_resistivity')
      if N == 1
         error('missing apparent resistivity conversion factor N');
      end
      b = N * b; % voltage -> apparent resistivity
   elseif strcmp(in, 'apparent_resistivity') && strcmp(out, 'voltage')
      if N == 1
         error('missing apparent resistivity conversion factor N');
      end
      b = N \ b; % apparent resistivity -> voltage
   else
      error('unknown conversion %s -> %s', in, out);
   end

function pass = do_unit_test
pass = 1;
% -------------
% ADAPTED FROM
% Create simulation data $Id: basic_iterative01.m 3829 2013-04-13 14:21:30Z bgrychtol $
%  http://eidors3d.sourceforge.net/tutorial/adv_image_reconst/basic_iterative.shtml
% 3D Model
imdl= mk_common_model('c2t4',16); % 576 elements
imdl.solve = 'inv_solve_abs_core';
imdl.reconst_type = 'absolute';
imdl.parameters.elem_working = 'log_conductivity';
imdl.inv_solve.calc_solution_error = 0;
%show_fem(imdl.fwd_model);
imgsrc= mk_image( imdl.fwd_model, 1);
% set homogeneous conductivity and simulate
vh=fwd_solve(imgsrc);
% set inhomogeneous conductivity and simulate
ctrs= interp_mesh(imdl.fwd_model);
x= ctrs(:,1); y= ctrs(:,2);
r1=sqrt((x+5).^2 + (y+5).^2); r2 = sqrt((x-50).^2 + (y-50).^2);
imgsrc.elem_data(r1<50)= 0.75;
imgsrc.elem_data(r2<30)= 1.8;
hh=figure; subplot(221); show_fem(imgsrc,1); axis tight; title('synthetic data');
% inhomogeneous data
vi=fwd_solve( imgsrc );
% add noise
%Add 30dB SNR noise to data
noise_level= std(vi.meas - vh.meas)/10^(30/20);
vi.meas = vi.meas + noise_level*randn(size(vi.meas));
% Reconstruct Images
img1= inv_solve(imdl, vi);
figure(hh); subplot(222); show_fem(img1,1); axis tight; title('#1 verbosity=default');
% -------------
disp('TEST: previous solved at default verbosity');
disp('TEST: now solve same at verbosity=0 --> should be silent');
imdl.parameters.verbose = 0;
imdl.parameters.meas_working = 'apparent_resistivity';
img2= inv_solve(imdl, vi);
figure(hh); subplot(223); show_fem(img2,1); axis tight; title('#2 verbosity=0');
if any(abs((img1.elem_data - img2.elem_data) / img1.elem_data) > eps*length(img1.elem_data)*1e3)
  disp('TEST:  img1 != img2 --> FAIL');
  pass = 0;
else
  disp('TEST:  img1 == img2 --> PASS');
end
% -------------
disp('TEST: try coarse2fine mapping');
imdl_tmp= mk_common_model('b2t4',16); % 256 elements
% convert fwd_model into rec_model
fmdl = imdl_tmp.fwd_model;
cmdl.type = fmdl.type;
cmdl.name = fmdl.name;
cmdl.nodes = fmdl.nodes;
cmdl.elems = fmdl.elems;
cmdl.gnd_node = 0;
% merge some elements
% TODO
% delete some other elements from around the boundary
ctrs= interp_mesh(cmdl);
x= ctrs(:,1); y= ctrs(:,2);
r=sqrt((x+5).^2 + (y+5).^2);
cmdl.elems(y-x > 150, :) = [];

% build c2f map
imdl.rec_model = cmdl;
c2f = mk_coarse_fine_mapping(imdl.fwd_model,cmdl);
imdl.fwd_model.coarse2fine = c2f;
% solve
%imdl.parameters.verbose = 10;
img3= inv_solve(imdl, vi);
%figure(hh); subplot(224); show_fem(cmdl,1); axis tight; title('#3 c2f');
figure(hh); subplot(223); show_fem(img3,1); axis tight; title('#3 c2f');
% check
e1 = c2f \ img1.elem_data;
e3 = img3.elem_data;
err = abs((e1 - e3) ./ e1);
err_thres = 0.40;
if any(err > err_thres) % maximum 15% error
  ni = find(err > err_thres);
  fprintf('TEST:  img1 != img3 --> FAIL max(err) = %0.2e on %d elements (thres=%0.2e)\n', ...
          max(err(ni)), length(ni), err_thres);
  pass = 0;
else
  disp('TEST:  img1 == img3 --> PASS');
end

imdl.parameters.elem_output = 'resistivity'; % resistivity output works
img4= inv_solve(imdl, vi);
figure(hh); subplot(224); show_fem(img4,1); axis tight; title('#4 c2f + resistivity out');
% check
e4 = 1./img4.elem_data;
err = abs((e1 - e4) ./ e1);
err_thres = 0.40;
if any(err > err_thres) % maximum 15% error
  ni = find(err > err_thres);
  fprintf('TEST:  img1 != img4 --> FAIL max(err) = %0.2e on %d elements (thres=%0.2e)\n', ...
          max(err(ni)), length(ni), err_thres);
  pass = 0;
else
  disp('TEST:  img1 == img4 --> PASS');
end


if pass
  disp('TEST: overall PASS');
else
  disp('TEST: overall FAIL');
end
