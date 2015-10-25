function img= inv_solve_abs_core( inv_model, data0);
%INV_SOLVE_ABS_CORE Absolute solver using a generic iterative algorithm
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data0      => EIT data
%
% This function is parameterized and uses function pointers where possible to
% allow its use as a general iterative solver framework. There are a large
% number of parameters and functions contained here. Sensible defaults are
% used throughout. You do not need to set every parameter.
%
% The solver operates as an absolute Gauss-Newton iterative solver by default.
% Wrapper functions are available to call this function in its various forms.
% Look forward to the "See also" section at the end of this help.
%
% Argument matrices to the internal functions (measurement inverse covariance,
% for example) are only calculated if required. Functions that are supplied to
% this INV_SOLVE_ABS_CORE must be able to survive being probed: they will have
% each parameter set to either 0 or 1 to determine if the function is sensitive
% to that argument. This works cleanly for most matrix multiplication based
% functions but for more abstract code, some handling of this behaviour may
% need to be implemented.
%
% In the following parameters, r_k is the current residual, r_{k-1} is the
% previous iteration's residual. k is the iteration count.
%
% Parameters denoted with a ** to the right of their default values are
% deprecated legacy parameters, some of which formerly existed under
% 'inv_model.parameters.*'.
%
% Parameters (inv_model.inv_solve_abs_core.*):
%   verbose (show progress)                (default 4)
%      0: quiet
%    >=1: print iteration count
%    >=2: print details as the algorithm progresses
%    >=3: plot residuals versus iteration count
%    >=4: plot result at each iteration, see show_fem
%    >=5: plot line search per iteration
%   plot_residuals                         (default 0)
%    plot residuals without verbose output
%   fig_prefix                       (default: <none>)
%    figure file prefix; figures not saved if <none>
%   fwd_solutions                          (default 0)
%    0: ignore
%    1: count fwd_solve(), generally the most
%       computationally expensive component of
%       the iterations
%   residual_func =             (default @GN_residual)
%    NOTE: @meas_residual exists to maintain
%    compatibility with some older code
%   max_iterations                        (default 10)  **
%   ntol (estimate of machine precision) (default eps)
%   tol (stop iter if r_k < tol)           (default 0)
%   dtol                              (default -0.01%)
%    stop iter if (r_k - r_{k-1}) < dtol AND
%                 k >= dtol_iter
%   dtol_iter                              (default 0)
%    apply dtol stopping criteria if k >= dtol_iter
%   min_value                           (default -inf)  **
%   max_value                           (default +inf)  **
%   line_optimize_func                (default <none>)  ** TODO
%     [next,fmin,res]=f(org,dx,data0,opt);
%     opt=line_optimize.* + objective_func
%     Deprecated, use line_search_func instead.
%   line_optimize.perturb
%                   (default line_search_args.perturb)  ** TODO
%     Deprecated, use line_search_args.perturb instead.
%   update_func                         (default TODO)  ** TODO
%     [img,opt]=f(org,next,dx,fmin,res,opt)
%     Deprecated, use <TODO> instead.
%   do_starting_estimate                   (default 1)  ** TODO
%     Deprecated, use <TODO> instead.
%   line_search_func       (default @line_search_onm2)
%   line_search_dv_func      (default @update_dv_core)
%   line_search_de_func      (default @update_de_core)
%   line_search_args.perturb
%                     (default [0 1/16 1/8 1/4 1/2 1])
%    line search for alpha by these steps along sx
%   line_search_args.plot                  (default 0)
%   c2f_background                         (default 0)
%    if > 0, this is additional elem_data
%    if a c2f map exists, the default is to decide
%    based on an estimate of c2f overlap whether a
%    background value is required. If a background is
%    required, it is added as the last element of that
%    type.
%   c2f_background_fixed                   (default 1)
%    hold the background estimate fixed or allow it
%    to vary as any other elem_data
%   elem_fixed                            (default [])
%    meas_select already handles selecting from the
%    valid measurements. we want the same for the
%    elem_data, so we only work on modifying the
%    legal values.
%    Note that c2f_background's elements are added to
%    this list if c2f_background_fixed == 1.
%   prior_data             (default to jacobian_bkgnd)
%    Sets the priors of type elem_prior. May be
%    scalar, per elem_prior, or match the working
%    length of each elem_data type. Note that for priors
%    using the c2f a background element may be added
%    to the end of that range when required; see
%    c2f_background.
%   elem_len                (default to all elem_data)
%    A cell array list of how many of each
%    elem_working there are in elem_data.
%      prior_data = { 32.1, 10*ones(10,1) };
%      elem_prior = {'conductivity', 'movement'};
%      elem_len = { 20001, 10 };
%   elem_prior                (default 'conductivity')
%    Input 'prior_data' type; immediately converted to
%    'elem_working' type before first iteration.
%   elem_working              (default 'conductivity')
%   elem_output               (default 'conductivity')
%    The working and output units for 'elem_data'.
%    Valid types are 'conductivity' and 'resistivity'
%    as plain units or with the prefix 'log_' or
%    'log10_'. Conversions are handled internally.
%    Scaling factors are applied to the Jacobian
%    (calculated in units of 'conductivity') as
%    appropriate.
%    If elem_working == elem_output, then no
%    conversions take place.
%    For multiple types, use cell array.
%    ex: elem_output = {'log_resistivity', 'movement'}
%   meas_input                     (default 'voltage')
%   meas_working                   (default 'voltage')
%    Similarly to elem_working/output, conversion
%    between 'voltage' and 'apparent_resistivity' and
%    their log/log10 varients are handled internally.
%    If meas_input == meas_working no conversions take
%    place. The normalization factor 'N' is calculated
%    if 'apparent_resistivity' is used.
%   update_img_func             (default: pass-through)
%    Called prior to calc_jacobian and update_dv.
%    Elements are converted to their "base types"
%    before this function is called. For example,
%    'log_resistivity' becomes 'conductivity'.
%    It is a hook to allow additional updates to the
%    model before the Jacobian, or a new set of
%    measurements are calculated via fwd_solve.
%   return_working_variables               (default: 0)
%    If 1, return the last working variables to the user
%     img.var.J   Jacobian
%     img.var.dx  descent direction
%     img.var.sx  search direction
%     img.var.alpha  line search result
%     img.var.beta   conjugation parameter
%     img.var.r   as:
%       [ residual, mesurement misfit, element misfit ]
%       with one row per iteration
%   show_fem                       (default: @show_fem)
%    Function with which to plot each iteration's
%    current parameters.
%
%   Signature for residual_func
%    [r,m,e] = f(dv, de, W, hp2RtR)
%   where
%    r   - the residual = m + e
%    m   - measurement error
%    e   - prior misfit
%    dv  - change in voltage
%    de  - change in image elements
%    W   - measurement inverse covarience matrix
%    hp2 - hyperparameter squared, see CALC_HYPERPARAMETER
%    RtR - regularization matrix squared --> hp2RtR = hp2*RtR
%
%   Signature for line_optimize_func
%    [alpha, img, dv, opt] = f(img, sx, data0, img0, N, W, hp2RtR, dv, opt)
%   where:
%    alpha - line search result
%    img   - the current image
%            (optional, recalculated if not available)
%    sx    - the search direction to which alpha should be applied
%    data0 - the true measurements     (dv = N*data - N*data0)
%    img0  - the image background (de = img - img0)
%    N     - a measurement normalization factor, N*dv
%    W     - measurement inverse covarience matrix
%    hp2   - hyperparameter squared, see CALC_HYPERPARAMETER
%    RtR   - regularization matrix squared --> hp2RtR = hp2*RtR
%    dv    - change in voltage
%            (optional, recalculated if not available)
%    opt   - additional arguments, updated at each call
%
%   Signature for line_search_dv_func
%    [dv, opt] = update_dv_core(img, data0, N, opt)
%   where:
%    dv    - change in voltage
%    opt   - additional arguments, updated at each call
%    data  - the estimated measurements
%    img   - the current image
%    data0 - the true measurements
%    N     - a measurement normalization factor, N*dv
%
%   Signature for line_search_de_func
%    de = f(img, img0, opt)
%   where:
%    de    - change in image elements
%    img   - the current image
%    img0  - the image background (de = img - img0)
%    opt   - additional arguments
%
%   Signature for calc_jacobian_scaling_func
%    S = f(x)
%   where:
%    S - to be used to scale the Jacobian
%    x - current img.elem_data in units of 'conductivity'
%
%   Signature for  update_img_func
%    img2 = f(img1, opt)
%   where
%    img1 - an input image, the current working image
%    img2 - a (potentially) modified version to be used
%    for the fwd_solve/Jacobian calculations
%
% NOTE that the default line search is very crude. For
% my test problems it seems to amount to an expensive grid
% search. Much more efficient line search algorithms exist
% and some fragments already are coded elsewhere in the
% EIDORS code-base.
%
% See also: INV_SOLVE_ABS_GN, INV_SOLVE_ABS_GN_LOGC,
%           INV_SOLVE_ABS_CG, INV_SOLVE_ABS_CG_LOGC,
%           LINE_SEARCH_O2, LINE_SEARCH_ONM2
%
% (C) 2010-2014 Alistair Boyle, Nolwenn Lesparre, Andy Adler, Bartłomiej Grychtol.
% License: GPL version 2 or version 3

% $Id$

%--------------------------
% UNIT_TEST?
if ischar(inv_model) && strcmp(inv_model,'UNIT_TEST') && (nargin == 1); do_unit_test; return; end
if ischar(inv_model) && strcmp(inv_model,'UNIT_TEST') && (nargin == 2); do_unit_test(data0); return; end

%--------------------------
opt = parse_options(inv_model);
if opt.verbose > 1
   fprintf('  verbose = %d\n', opt.verbose);
end
%if opt.do_starting_estimate
%    img = initial_estimate( inv_model, data0 ); % TODO
%%%    AB->NL this is Nolwenn's homogeneous estimate...
%%%    calc_background_resistivity is my version of this code
%%%    that is working for my data set
%else
[inv_model, opt] = append_c2f_background(inv_model, opt);
% calc_jacobian_bkgnd, used by mk_image does not understand
% the course-to-fine mapping and explodes when it is fed a
% prior based on the coarse model. Here we give that
% function something it can swallow, then create then plug
% in the correct prior afterwards.
if isfield(inv_model, 'jacobian_bkgnd')
  inv_model = rmfield(inv_model,'jacobian_bkgnd');
end
inv_model.jacobian_bkgnd.value = 1;
img = mk_image( inv_model );
img.inv_model = inv_model; % stash the inverse model
img = data_mapper(img); % move data from whatever 'params' to img.elem_data

% insert the prior data
img = init_elem_data(img, opt);

% map data and measurements to working types
%  convert elem_data
img = map_img(img, opt.elem_working);
%  convert measurement data
if ~isstruct(data0)
   d = data0;
   data0 = struct;
   data0.meas = d;
   data0.type = 'data';
end
data0.current_params = opt.meas_input;

% precalculate some of our matrices if required
W  = init_meas_icov(inv_model, opt);
[N, dN] = init_normalization(inv_model.fwd_model, data0, opt);

% now get on with
img0 = img;
hp2RtR = 0; alpha = 0; k = 0; sx = 0; r = 0; stop = 0; % general init
residuals = zeros(opt.max_iterations,3); % for residuals plots
dxp = 0; % previous step's slope was... nothing
[dv, opt] = update_dv([], img, data0, N, opt);
de = update_de([], img, img0, opt);
if opt.verbose >= 5 % we only save the measurements at each iteration if we are being verbose
  dvall = ones(size(data0.meas,1),opt.max_iterations+1)*NaN;
end
while 1
  if opt.verbose > 1
     if k == 0
        fprintf('  iteration start up\n')
     else
        fprintf('  iteration %d\n', k)
     end
  end

  % calculate the Jacobian
  %  - Jacobian before RtR because it is needed for Noser prior
  [J, opt] = update_jacobian(img, dN, k, opt);

  % update RtR, if required (depends on prior)
  hp2RtR = update_hp2RtR(inv_model, J, k, img, opt);

  % determine the next search direction sx
  %  dx is specific to the algorithm, generally "downhill"
  dx = update_dx(J, W, hp2RtR, dv, de, opt);
  % choose beta, beta=0 unless doing Conjugate Gradient
  beta = update_beta(dx, dxp, sx, opt);
  beta_all(k+1)=beta; % save for debug
  % sx_k = dx_k + beta * sx_{k-1}
  sx = update_sx(dx, beta, sx, opt);
  if k ~= 0
     dxp = dx; % saved for next iteration if using beta
  end

  % plot SVD of Jacobian (before and after regularization)
  plot_svd_elem(J, W, hp2RtR, k, sx, dx, img, opt);

  % line search for alpha, leaving the final selection as img
  % x_n = img.elem_data
  % x_{n+1} = x_n + \alpha sx
  % img.elem_data = x_{n+1}
  [alpha, img, dv, opt] = update_alpha(img, sx, data0, img0, N, W, hp2RtR, k, dv, opt);
  alpha_all(k+1) = alpha;
  % fix max/min values for x, clears dx if limits are hit, where
  % a cleared dv will trigger a recalculation of dv at the next update_dv()
  [img, dv] = update_img_using_limits(img, img0, data0, N, dv, opt);

  % update change in element data from the prior de and
  % the measurement error dv
  [dv, opt] = update_dv(dv, img, data0, N, opt);
  de = update_de(de, img, img0, opt);
  if opt.verbose >= 5
    dvall(:,k+1) = dv;
    show_meas_err(dvall, data0, k+1, N, W, opt);
  end
  show_fem_iter(k, img, inv_model, stop, opt);

  % now find the residual, quit if we're done
  [stop, k, r, img] = update_residual(dv, img, de, W, hp2RtR, k, r, alpha, sx, opt);
  if stop
     if stop == -1
        alpha_all(k) = 0;
     end
     break;
  end
end
[img, opt] = strip_c2f_background(img, opt);
% check we're returning the right size of data
if isfield(inv_model, 'rec_model')
  img.fwd_model = inv_model.rec_model;
end
if opt.verbose > 1
   if k==1; itrs=''; else itrs='s'; end
   fprintf('  %d fwd_solves required for this solution in %d iteration%s\n', ...
           opt.fwd_solutions, k, itrs);
end
% convert data for output
img = map_img(img, opt.elem_output);
img.meas_err = dv;
if opt.return_working_variables
  img.inv_solve_abs_core.J = J;
  img.inv_solve_abs_core.dx = dx;
  img.inv_solve_abs_core.sx = sx;
  img.inv_solve_abs_core.alpha = alpha_all;
  img.inv_solve_abs_core.beta = beta_all;
  img.inv_solve_abs_core.k = k;
  img.inv_solve_abs_core.r = r(1:(k+1),:); % trim r to n-iterations' rows
  img.inv_solve_abs_core.N = N;
  img.inv_solve_abs_core.W = W;
  img.inv_solve_abs_core.hp2RtR = hp2RtR;
  img.inv_solve_abs_core.dv = dv;
  img.inv_solve_abs_core.de = de;
  if opt.verbose >= 5
    img.inv_solve_abs_core.dvall = dvall;
  end
end
%img = data_mapper(img, 1); % move data from img.elem_data to whatever 'params'

function show_meas_err(dvall, data0, k, N, W, opt)
   clf;
   subplot(211); bar(dvall(:,k)); ylabel(sprintf('dv_k [%s]',opt.meas_working)); xlabel('meas #'); title(sprintf('iter %d',k));
   subplot(212); bar(map_meas(dvall(:,k),N,opt.meas_working, 'voltage')); ylabel('dv_k [V]'); xlabel('meas #'); title('');
   drawnow;
   if isfield(opt,'fig_prefix')
      print('-dpdf',sprintf('%s-meas_err%d',opt.fig_prefix,k));
      print('-dpng',sprintf('%s-meas_err%d',opt.fig_prefix,k));
      saveas(gcf,sprintf('%s-meas_err%d.fig',opt.fig_prefix,k));
   end
   drawnow;

function img = init_elem_data(img, opt)
  if opt.verbose > 1
    fprintf('  setting prior elem_data\n');
  end
  ne2 = 0; % init
  img.elem_data = zeros(sum([opt.elem_len{:}]),1); % preallocate
  for i=1:length(opt.elem_prior)
    ne1 = ne2+1; % next start idx ne1
    ne2 = ne1+opt.elem_len{i}-1; % this set ends at idx ne2
    if opt.verbose > 1
      if length(opt.prior_data{i}) == 1
        fprintf('    %d x %s: %0.1f\n',opt.elem_len{i},opt.elem_prior{i}, opt.prior_data{i});
      else
        fprintf('    %d x %s: ...\n',opt.elem_len{i},opt.elem_prior{i});
        if length(opt.prior_data{i}) ~= opt.elem_len{i}
           error(sprintf('expected %d elem, got %d elem in elem_prior', ...
                         opt.elem_len{i}, length(opt.prior_data{i})));
        end
      end
    end
    img.params_sel(i) = {ne1:ne2};
    img.elem_data(img.params_sel{i}) = opt.prior_data{i};
  end
  img.current_params = opt.elem_prior;

function W = init_meas_icov(inv_model, opt)
   W = 1;
   if opt.calc_meas_icov
      if opt.verbose > 1
         disp('  calc measurement inverse covariance W');
      end
      W   = calc_meas_icov( inv_model );
   end
   err_if_inf_or_nan(W, 'init_meas_icov');

function [N, dN] = init_normalization(fmdl, data0, opt)
   % precalculate the normalization of the data if required (apparent resistivity)
   N = 1;
   dN = 1;
   vh1.meas = 1;
   if ~ischar(opt.meas_input) || ~ischar(opt.meas_working)
      error('expected strings for meas_input and meas_working');
   end
   go =       any(strcmp({opt.meas_input, opt.meas_working},'apparent_resistivity'));
   go = go || any(strcmp({opt.meas_input, opt.meas_working},'log_apparent_resistivity'));
   go = go || any(strcmp({opt.meas_input, opt.meas_working},'log10_apparent_resistivity'));
   if go
      if opt.verbose > 1
         disp(['  calc measurement normalization matrix N (voltage -> ' opt.meas_working ')']);
      end
      % calculate geometric factor for apparent_resitivity conversions
      img1 = mk_image(fmdl,1);
      vh1  = fwd_solve(img1);
      N    = spdiag(1./vh1.meas);
      err_if_inf_or_nan(N,  'init_normalization: N');
   end
   if go && (opt.verbose > 1)
      disp(['  calc Jacobian normalization matrix   dN (voltage -> ' opt.meas_working ')']);
   end
   % calculate the normalization factor for the Jacobian
   data0 = map_meas_struct(data0, N, 'voltage'); % to voltage
   switch opt.meas_working
      case 'apparent_resistivity'
         dN = da_dv(data0.meas, vh1.meas);
      case 'log_apparent_resistivity'
         dN = dloga_dv(data0.meas, vh1.meas);
      case 'log10_apparent_resistivity'
         dN = dlog10a_dv(data0.meas, vh1.meas);
      case 'voltage'
         dN = dv_dv(data0.meas, vh1.meas);
      case 'log_voltage'
         dN = dlogv_dv(data0.meas, vh1.meas);
      case 'log10_voltage'
         dN = dlog10v_dv(data0.meas, vh1.meas);
      otherwise
         error('hmm');
   end
   err_if_inf_or_nan(dN, 'init_normalization: dN');

% r_km1: previous residual, if its the first iteration r_km1 = inf
% r_k: new residual
function [stop, k, r, img] = update_residual(dv, img, de, W, hp2RtR, k, r, alpha, sx, opt)
  stop = 0;

  % update residual estimate
  if k == 0
     r = ones(opt.max_iterations, 3)*NaN;
     r_km1 = inf;
  else
     r_km1 = r(k, 1);
  end
  [r_k m_k e_k] = feval(opt.residual_func, dv, de, W, hp2RtR);
  % save residual for next iteration
  r(k+1,1:3) = [r_k m_k e_k];

  % now do something with that information
  if opt.verbose > 1
     if k == 0
        fprintf('    stop @ max iter = %d, tol = %0.3g (%0.3g%%), dtol = %0.3g%% (after %d iter)\n', ...
                opt.max_iterations, opt.tol, opt.tol/r_k*100, opt.dtol*100, opt.dtol_iter);
        fprintf('    calc residual, r=%0.3g\n', r_k);
     else
        fprintf('    calc residual\n');
        fprintf('      r =%0.3g (%0.03g%%)\n', r_k, r_k/r(1)*100);
        dr = (r_k - r_km1);
        fprintf('      dr=%0.3g (%0.3g%%)\n', dr, dr/r_km1*100);
     end
  end
  if opt.plot_residuals
     %         optimization_criteria, data misfit, roughness
     r(k+1,2:3) = [(dv'*dv)/2 (de'*de)/2];
     if k > 0
        clf;
        x = 1:(k+1);
        y = r(x, :);
        y = y ./ repmat(max(y,[],1),size(y,1),1) * 100;
        plot(x-1, y, 'o-', 'linewidth', 2, 'MarkerSize', 10);
        title('residuals');
        axis tight;
        ylabel('residual (% of max)');
        xlabel('iteration');
        set(gca, 'xtick', x);
        set(gca, 'xlim', [0 max(x)-1]);
        legend('residual','meas. misfit','prior misfit');
        legend('Location', 'EastOutside');
        drawnow;
        if isfield(opt,'fig_prefix')
           print('-dpdf',sprintf('%s-r%d',opt.fig_prefix,k));
           print('-dpng',sprintf('%s-r%d',opt.fig_prefix,k));
           saveas(gcf,sprintf('%s-r%d.fig',opt.fig_prefix,k));
        end
     end
  end

  % evaluate stopping criteria
  if r_k > r_km1 % bad step
     if opt.verbose > 1
        fprintf('  terminated at iteration %d (bad step, returning previous iteration''s result)\n',k);
     end
     img.elem_data = img.elem_data - alpha * sx; % undo the last step
     stop = -1;
  elseif k >= opt.max_iterations
     if opt.verbose > 1
        fprintf('  terminated at iteration %d (max iterations)\n',k);
     end
     stop = 1;
  elseif r_k < opt.tol + opt.ntol
     if opt.verbose > 1
        fprintf('  terminated at iteration %d\n',k);
        fprintf('    residual tolerance (%0.3g) achieved\n', opt.tol + opt.ntol);
     end
     stop = 1;
  elseif (k >= opt.dtol_iter) && ((r_k - r_km1)/r_km1 > opt.dtol + 2*opt.ntol)
     if opt.verbose > 1
        fprintf('  terminated at iteration %d (iterations not improving)\n', k);
        fprintf('    residual slope tolerance (%0.3g%%) exceeded\n', (opt.dtol + 2*opt.ntol)*100);
     end
     stop = 1;
  end
  if ~stop
     % update iteration count
     k = k+1;
  end

% for Conjugate Gradient, else beta = 0
%  dx_k, dx_{k-1}, sx_{k-1}
function beta = update_beta(dx_k, dx_km1, sx_km1, opt);
   if isfield(opt, 'beta_func')
      if opt.verbose > 1
         try beta_str = func2str(opt.beta_func);
         catch
            try beta_str = opt.beta_func;
            catch beta_str = 'unknown';
            end
         end
      end
      beta= feval(opt.beta_func, dx_k, dx_km1, sx_km1);
   else
     beta_str = '<none>';
     beta = 0;
   end
   if opt.verbose > 1
      str = sprintf('    calc beta (%s)=%0.3f\n', beta_str, beta);
   end

% update the search direction
% for Gauss-Newton
%   sx_k = dx_k
% for Conjugate-Gradient
%   sx_k = dx_k + beta * sx_{k-1}
function sx = update_sx(dx, beta, sx_km1, opt);
   sx = dx + beta * sx_km1;
   if(opt.verbose > 1)
      nsx = norm(sx);
      nsxk = norm(sx_km1);
      fprintf( '    update step dx, beta=%0.3g, ||dx||=%0.3g\n', beta, nsx);
      if nsxk ~= 0
         fprintf( '      acceleration     d||dx||=%0.3g\n', nsx-nsxk);
         % ||ddx|| = chord_len = 2 sin(theta/2)
         ddx = norm(sx/nsx-sx_km1/nsxk);
         fprintf( '      direction change ||ddx||=%0.3g (%0.3g°)\n', ddx, 2*asind(ddx/2)); 
      end
   end

% this function constructs the blockwise RtR matrix
function hp2RtR = update_hp2RtR(inv_model, J, k, img, opt)
   if k==0 % first the start up iteration use the initial hyperparameter
      k=1;
   end
   % TODO sometimes (with Noser?) this requires the Jacobian, could this be done more efficiently?
   % add a test function to determine if img.elem_data affects RtR, skip this if independant
   % TODO we could detect in the opt_parsing whether the calc_RtR_prior depends on 'x' and skip this if no
   if ~opt.calc_RtR_prior
      error('no RtR calculation mechanism, set imdl.inv_solve_abs_core.RtR_prior or imdl.RtR_prior');
   end
   if opt.verbose > 1
      disp('    calc hp^2 R^t R');
   end
   hp2  = calc_hyperparameter( inv_model ); % = \lambda^2
   net = sum([opt.elem_len{:}]); % Number of Elements, Total
   RtR = sparse(net,net); % init RtR = sparse(zeros(net,net));
   esi = 0; eei = 0; % element start, element end
   for i = 1:size(opt.RtR_prior,1) % row i
      esi = eei + 1;
      eei = eei + opt.elem_len{i};
      esj = 0; eej = 0; % element start, element end
      for j = 1:size(opt.RtR_prior,2) % column j
         esj = eej + 1;
         eej = eej + opt.elem_len{j};
         if isempty(opt.RtR_prior{i,j}) % null entries
            continue; % no need to explicitly create zero block matrices
         end

         % select a hyperparameter, potentially, per iteration
         % if we're at the end of the list, select the last entry
         hp=opt.hyperparameter{i,j};
         if length(hp) > k
            hp=hp(k);
         else
            hp=hp(end);
         end

         if opt.verbose > 1
            try RtR_str = func2str(opt.RtR_prior{i,j});
            catch
               try RtR_str = opt.RtR_prior{i,j};
               catch RtR_str = 'unknown';
               end
            end
            fprintf('      {%d,%d} regularization RtR (%s), ne=%dx%d, hp=%0.4g\n', i,j,RtR_str,eei-esi+1,eej-esj+1,hp*sqrt(hp2));
         end
         imgt = map_img(img, opt.elem_working{i});
         inv_modelt = inv_model;
         inv_modelt.RtR_prior = opt.RtR_prior{i,j};
         RtR(esi:eei, esj:eej) = hp.^2 * calc_RtR_prior_wrapper(inv_modelt, imgt, opt);
      end
   end
   hp2RtR = hp2*RtR;

function plot_svd_elem(J, W, hp2RtR, k, sx, dx, img, opt)
   if(opt.verbose >= 5)
      % and try a show_fem with the pixel search direction
      clf;
      imgb=img;
      imgb.elem_data = dx;
      imgb.current_params = opt.elem_working;
      if isfield(imgb.inv_model,'rec_model')
         imgb.fwd_model = imgb.inv_model.rec_model;
      end
      feval(opt.show_fem,imgb,[]);
      title(sprintf('dx @ iter=%d',k));
      drawnow;
      if isfield(opt,'fig_prefix')
         print('-dpng',sprintf('%s-dx%d',opt.fig_prefix,k));
         print('-dpdf',sprintf('%s-dx%d',opt.fig_prefix,k));
         saveas(gcf,sprintf('%s-dx%d.fig',opt.fig_prefix,k));
      end
   end
   if opt.verbose < 8
      return; % do nothing if not verbose
   end
   % canonicalize the structures so we don't have to deal with a bunch of scenarios below
   if ~isfield(img, 'params_sel')
      img.params_sel = {1:length(img.elem_data)};
   end
   if ~isfield(img, 'current_params')
      img.current_params = 'conductivity';
   end
   if ~iscell(img.current_params)
      img.current_params = {img.current_params};
   end
   % go
   cols=length(opt.elem_working);
   if norm(sx - dx) < range(dx)/max(dx)*0.01 % sx and dx are within 1%
      rows=2;
   else
      rows=3;
   end
   clf; % individual SVD plots
   for i=1:cols
      if 1 % if Tikhonov
         hp=opt.hyperparameter{i};
         if k ~= 0 && k < length(hp)
            hp = hp(k);
         else
            hp = hp(end);
         end
      else
         hp = [];
      end
      sel=img.params_sel{i};
      str=strrep(img.current_params{i},'_',' ');
      plot_svd(J(:,sel), W, hp2RtR(sel,sel), k, hp); xlabel(str);
      drawnow;
      if isfield(opt,'fig_prefix')
         print('-dpdf',sprintf('%s-svd%d-%s',opt.fig_prefix,k,img.current_params{i}));
         print('-dpng',sprintf('%s-svd%d-%s',opt.fig_prefix,k,img.current_params{i}));
         saveas(gcf,sprintf('%s-svd%d-%s.fig',opt.fig_prefix,k,img.current_params{i}));
      end
   end
   clf; % combo plot
   for i=1:cols
      if 1 % if Tikhonov
         hp=opt.hyperparameter{i};
         if k ~= 0 && k < length(hp)
            hp = hp(k);
         else
            hp = hp(end);
         end
      else
         hp = [];
      end
      subplot(rows,cols,i);
      sel=img.params_sel{i};
      str=strrep(img.current_params{i},'_',' ');
      plot_svd(J(:,sel), W, hp2RtR(sel,sel), k, hp); xlabel(str);
      subplot(rows,cols,cols+i);
      bar(dx(sel)); ylabel(['dx: ' str]);
      if rows > 2
         subplot(rows,cols,2*cols+i);
         bar(sx(sel)); ylabel(['sx: ' str]);
      end
   end
   drawnow;
   if isfield(opt,'fig_prefix')
      print('-dpdf',sprintf('%s-svd%d',opt.fig_prefix,k));
      print('-dpng',sprintf('%s-svd%d',opt.fig_prefix,k));
      saveas(gcf,sprintf('%s-svd%d.fig',opt.fig_prefix,k));
   end

function plot_svd(J, W, hp2RtR, k, hp)
   if nargin < 5
      hp = [];
   end
   % calculate the singular values before and after regularization
   [~,s1,~]=svd(J'*W*J); s1=sqrt(diag(s1));
   [~,s2,~]=svd(J'*W*J + hp2RtR); s2=sqrt(diag(s2));
   h=semilogy(s1,'bx'); axis tight; set(h,'LineWidth',2);
   hold on; h=semilogy(s2,'go'); axis tight; set(h,'LineWidth',2); hold off;
   xlabel('k'); ylabel('value \sigma');
   title(sprintf('singular values of J at iteration %d',k));
   legend('J^T J', 'J^T J + \lambda^2 R^T R'); legend location best;
   % line for \lambda
%     if regularization == 2 % Noser
%        hp_scaled = hp*sqrt(norm(full(RtR)));
%        h=line([1 length(s1)],[hp_scaled hp_scaled]);
%        text(length(s1)/2,hp_scaled*0.9,sprintf('\\lambda ||R^T R||^{%0.1f}= %0.4g; \\lambda = %0.4g', noser_p, hp_scaled, hp));
%        fprintf('  affecting %d of %d singular values\k', length(find(s1<hp_scaled)), length(s1));
%     else % Tikhonov
   if length(hp)==1
        h=line([1 length(s1)],[hp hp]);
        ly=10^(log10(hp)-0.05*range(log10([s1;s2])));
        text(length(s1)/2,ly,sprintf('\\lambda = %0.4g', hp));
%       fprintf('  affecting %d of %d singular values\k', length(find(s1<hp)), length(s1));
     end
     set(h,'LineStyle','-.'); set(h,'LineWidth',2);
   set(gca,'YMinorTick','on', 'YMinorGrid', 'on', 'YGrid', 'on');


% TODO this function is one giant HACK around broken RtR generation with c2f matrices
function RtR = calc_RtR_prior_wrapper(inv_model, img, opt)
   RtR = calc_RtR_prior( inv_model );
   if size(RtR,1) < length(img.elem_data)
     ne = length(img.elem_data) - size(RtR,1);
     % we are correcting for the added background element
     for i=1:ne
       RtR(end+1:end+1, end+1:end+1) = RtR(1,1);
     end
     if opt.verbose > 1
        fprintf('      c2f: adjusting RtR by appending %d rows/cols\n', ne);
        disp(   '      TODO move this fix, or something like it to calc_RtR_prior -- this fix is a quick HACK to get things to run...');
     end
   end

% opt is only updated for the fwd_solve count
function [J, opt] = update_jacobian(img, dN, k, opt)
   k=k+1;
   base_types = map_img_base_types(img);
   imgb = map_img(img, base_types);
   imgb = feval(opt.update_img_func, imgb, opt);
   % if the electrodes/geometry moved, we need to recalculate dN if it depends on vh
   % note that only apparent_resisitivity needs vh; all others depend on data0 measurements
   if any(strcmp(map_img_base_types(img), 'movement')) && any(strcmp(opt.meas_working, 'apparent_resistivity'))
      imgh = map_img(imgb, 'conductivity'); % drop everything but conductivity
      imgh.elem_data = imgh.elem_data*0 +1; % conductivity = 1
      vh = fwd_solve(imgh); vh = vh.meas;
      dN = da_dv(1,vh); % = diag(1/vh)
      opt.fwd_solutions = opt.fwd_solutions +1;
   end
   ee = 0; % element select, init
   pp = fwd_model_parameters(imgb.fwd_model);
   J = zeros(pp.n_meas,sum([opt.elem_len{:}]));
   for i=1:length(opt.jacobian)
     if(opt.verbose > 1)
        try J_str = func2str(opt.jacobian{i});
        catch J_str = opt.jacobian{i};
        end
        if i == 1 fprintf('    calc Jacobian J(x) = ');
        else      fprintf('                       + '); end
        fprintf('(%s,', J_str);
     end
     % start and end of these Jacobian columns
     es = ee+1;
     ee = es+opt.elem_len{i}-1;
     % scaling if we are working in something other than direct conductivity
     S = feval(opt.calc_jacobian_scaling_func{i}, imgb.elem_data(es:ee)); % chain rule
     % finalize the jacobian
     % Note that if a normalization (i.e. apparent_resistivity) has been applied
     % to the measurements, it needs to be applied to the Jacobian as well!
     imgt = imgb;
     if  strcmp(base_types{i}, 'conductivity') % make legacy jacobian calculators happy... only conductivity on imgt.elem_data
        imgt = map_img(img, 'conductivity');
     end
     imgt.fwd_model.jacobian = opt.jacobian{i};
     Jn = calc_jacobian( imgt ); % unscaled natural units (i.e. conductivity)
     J(:,es:ee) = dN * Jn * S; % scaled and normalized
     if opt.verbose > 1
        tmp = zeros(1,size(J,2));
        tmp(es:ee) = 1;
        tmp(opt.elem_fixed) = 0;
        fprintf(' %d DoF, %d meas, %s)\n', sum(tmp), size(J,1), func2str(opt.calc_jacobian_scaling_func{i}));
     end
     if opt.verbose >= 5
        clf;
        t=axes('Position',[0 0 1 1],'Visible','off'); % something to put our title on after we're done
        text(0.03,0.1,sprintf('update\\_jacobian (%s), iter=%d', strrep(J_str,'_','\_'), k),'FontSize',20,'Rotation',90);
        for y=0:1
           if y == 0; D = Jn; else D = J(:,es:ee); end
           axes('units', 'normalized', 'position', [ 0.13 0.62-y/2 0.8 0.3 ]);
           imagesc(D);
           if y == 0; ylabel('meas (1)'); xlabel(['elem (' strrep(base_types{i},'_','\_') ')']);
           else       ylabel('meas (dN)'); xlabel(['elem (' strrep(opt.elem_working{i},'_','\_') ')']);
           end
           os = get(gca, 'Position'); c=colorbar('southoutside'); % colorbar start...
           set(gca, 'Position', os); % fix STUPID colorbar resizing
           % reduce height, this has to be done after the axes fix or STUPID matlab messes things up real good
           cP = get(c,'Position'); set(c,'Position', [0.13    0.54-y/2    0.8    0.010]);
           axes('units', 'normalized', 'position', [ 0.93 0.62-y/2 0.05 0.3 ]);
           barh(sqrt(sum(D.^2,2))); axis tight; axis ij; set(gca, 'ytick', [], 'yticklabel', []);
           axes('units', 'normalized', 'position', [ 0.13 0.92-y/2 0.8 0.05 ]);
           bar(sqrt(sum(D.^2,1))); axis tight; set(gca, 'xtick', [], 'xticklabel', []);
        end
        drawnow;
        if isfield(opt,'fig_prefix')
           print('-dpng',sprintf('%s-J%d-%s',opt.fig_prefix,k,strrep(J_str,'_','')));
           print('-dpdf',sprintf('%s-J%d-%s',opt.fig_prefix,k,strrep(J_str,'_','')));
           saveas(gcf,sprintf('%s-J%d-%s.fig',opt.fig_prefix,k,strrep(J_str,'_','')));
        end
        % and try a show_fem with the pixel sensitivity
        clf;
        imgb.elem_data = log(sqrt(sum(D.^2,1)));
        imgb.current_params = [ opt.elem_working{:} ' sensitivity'];
        if isfield(imgb.inv_model,'rec_model')
           imgb.fwd_model = imgb.inv_model.rec_model;
        end
        feval(opt.show_fem,imgb,[]);
        title(sprintf('sensitivity @ iter=%d',k));
        drawnow;
        if isfield(opt,'fig_prefix')
           print('-dpng',sprintf('%s-Js%d-%s',opt.fig_prefix,k,strrep(J_str,'_','')));
           print('-dpdf',sprintf('%s-Js%d-%s',opt.fig_prefix,k,strrep(J_str,'_','')));
           saveas(gcf,sprintf('%s-Js%d-%s.fig',opt.fig_prefix,k,strrep(J_str,'_','')));
        end
     end
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
   S = spdiag(x);
function S = dx_dlog10x(x);
   S = spdiag(x * log(10));
% resistivity 'y'
% d x     d x   -1
% ----- = --- = ---, y = 1/x --> -(x^2)
% d 1/x   d y   y^2
function S = dx_dy(x);
   S = spdiag(-(x.^2));
% then build the log versions of conductivity by combining chain rule products
function S = dx_dlogy(x);
%   S = dx_dy(x) * dy_dlogy(x);
%     = -(x^2) * 1/x = -x
   S = spdiag(-x);
function S = dx_dlog10y(x);
%   S = dx_dy(x) * dy_dlog10y(x);
%     = -(x^2) * 1/(ln(10) x) = -x / ln(10)
   S = spdiag(-x/log(10));
% ... some renaming to make things understandable above: x = 1/y
%function S = dy_dlogy(x);
%   S = dx_dlogx(1./x);
%function S = dy_dlog10y(x);
%   S = dx_dlog10x(1./x);
% -------------------------------------------------
% apparent_resistivity 'a' versus voltage 'x'
% d a    1  d v    1         v
% --- = --- --- = --- ; a = ---
% d v    vh d v    vh        vh
% log_apparent_resistivity
% d loga   d loga d a    1   1     vh  1     1
% ------ = ------ --- = --- --- = --- --- = ---
% d v       d a   d v    a   vh    v   vh    v
function dN = da_dv(v,vh)
   dN = spdiag(1./vh); % N == dN for apparent_resistivity
function dN = dloga_dv(v,vh)
   dN = spdiag(1./v);
function dN = dlog10a_dv(v,vh)
   dN = spdiag( 1./(v * log(10)) );
function dN = dv_dv(v,vh)
   dN = 1;
function dN = dlogv_dv(v,vh) % same as dloga_dv
   dN = dloga_dv(v,vh);
function dN = dlog10v_dv(v,vh) % same as dlog10a_dv
   dN = dlog10a_dv(v, vh);
% -------------------------------------------------


function [alpha, img, dv, opt] = update_alpha(img, sx, data0, img0, N, W, hp2RtR, k, dv, opt)
  if k == 0 % first iteration, just setting up, no line search happens
     alpha = 0;
     return;
  end

  if(opt.verbose > 1)
     try ls_str = func2str(opt.line_search_func);
     catch ls_str = opt.line_search_func;
     end
     fprintf('    line search, alpha = %s\n', ls_str);
  end

  % some sanity checks before we feed this information to the line search
  err_if_inf_or_nan(sx, 'sx (pre-line search)');
  err_if_inf_or_nan(img.elem_data, 'img.elem_data (pre-line search)');

  if any(size(img.elem_data) ~= size(sx))
     error(sprintf('mismatch on elem_data[%d,%d] vs. sx[%d,%d] vector sizes, check c2f_background_fixed',size(img.elem_data), size(sx)));
  end

  [alpha, imgo, dv, opto] = feval(opt.line_search_func, img, sx, data0, img0, N, W, hp2RtR, dv, opt);
  if ~isempty(imgo)
     img = imgo;
  else
     img.elem_data = img.elem_data + alpha*sx;
  end
  if ~isempty(opto)
     opt = opto;
  end

  if(opt.verbose > 1)
     fprintf('      selected alpha=%0.3g\n', alpha);
  end

  if (alpha == 0) && (k == 1)
    error('first iteration failed to advance solution');
  end

function err_if_inf_or_nan(x, str);
  if any(any(isnan(x) | isinf(x)))
      error(sprintf('bad %s (%d NaN, %d Inf of %d)', ...
                    str, ...
                    length(find(isnan(x))), ...
                    length(find(isinf(x))), ...
                    length(x(:))));
  end


function [img, dv] = update_img_using_limits(img, img0, data0, N, dv, opt)
  % fix max/min values for x
  if opt.max_value ~= +inf
     lih = find(img.elem_data > opt.max_value);
     img.elem_data(lih) = opt.max_value;
     if opt.verbose > 1
        fprintf('    limit max(x)=%g for %d elements\n', opt.max_value, length(lih));
     end
     dv = []; % dv is now invalid since we changed the conductivity
  end
  if opt.min_value ~= -inf
     lil = find(img.elem_data < opt.min_value);
     img.elem_data(lil) = opt.min_value;
     if opt.verbose > 1
        fprintf('    limit min(x)=%g for %d elements\n', opt.min_value, length(lil));
     end
     dv = []; % dv is now invalid since we changed the conductivity
  end
  % update voltage change estimate if the limit operation changed the img data
  [dv, opt] = update_dv(dv, img, data0, N, opt, '(dv out-of-date)');

function  de = update_de(de, img, img0, opt)
   img0 = map_img(img0, opt.elem_working);
   img  = map_img(img,  opt.elem_working);
   err_if_inf_or_nan(img0.elem_data, 'de img0');
   err_if_inf_or_nan(img.elem_data,  'de img');
   % probably not the most robust check for whether this is the first update
   % but this ensures that we get exactly zero for the first iteration and not
   % a set of values that has numeric floating point errors that are nearly zero
   if isempty(de) % first iteration
      % data hasn't changed yet!
      de = zeros(size(img0.elem_data));
   else
      de = img0.elem_data - img.elem_data;
   end
   de(opt.elem_fixed) = 0; % TODO is this redundant... delete me?
   err_if_inf_or_nan(de, 'de out');

function [dv, opt] = update_dv(dv, img, data0, N, opt, reason)
   % estimate current error as a residual
   if ~isempty(dv) % need to calculate dv...
      return;
   end
   if nargin < 7
      reason = '';
   end
   if opt.verbose > 1
      disp(['    fwd_solve b=Ax ', reason]);
   end
   [dv, opt, err] = update_dv_core(img, data0, N, opt);
% TODO AB inject the img.error here, so it doesn't need to be recalculated when calc_solution_error=1
%   img.error = err;

function data = map_meas_struct(data, N, out)
   try   current_meas_params = data.current_params;
   catch current_meas_params = 'voltage';
   end
   data.meas = map_meas(data.meas, N, current_meas_params, out);
   data.current_params = out;
   err_if_inf_or_nan(data.meas, 'dv meas');

% also used by the line search as opt.line_search_dv_func
function [dv, opt, err] = update_dv_core(img, data0, N, opt)
   data0 = map_meas_struct(data0, N, 'voltage');
   img = map_img(img, map_img_base_types(img));
   img = feval(opt.update_img_func, img, opt);
   img = map_img(img, 'conductivity'); % drop everything but conductivity
   % if the electrodes/geometry moved, we need to recalculate N if it's being used
   if any(any(N ~= 1)) && any(strcmp(map_img_base_types(img), 'movement'))
      % note: data0 is mapped back to 'voltage' before N is modified
      imgh=img; imgh.elem_data = imgh.elem_data*0 +1; % conductivity = 1
      vh = fwd_solve(imgh); vh = vh.meas;
      N = spdiag(1./vh);
      opt.fwd_solutions = opt.fwd_solutions +1;
   end
   data = fwd_solve(img);
   opt.fwd_solutions = opt.fwd_solutions +1;
   dv = calc_difference_data(data, data0, img.fwd_model);
   if nargout >= 3
      err = norm(dv)/norm(data0.meas);
   else
      err = NaN;
   end
   dv = map_meas(dv, N, 'voltage', opt.meas_working);
   err_if_inf_or_nan(dv, 'dv out');

function show_fem_iter(k, img, inv_model, stop, opt)
  if (opt.verbose < 4) || (stop == -1)
     return; % if verbosity is low OR we're dropping the last iteration because it was bad, do nothing
  end
  if opt.verbose > 1
     str=opt.show_fem;
     if isa(str,'function_handle')
        str=func2str(str);
     end
     disp(['    ' str '()']);
  end
  if isequal(opt.show_fem,@show_fem) % opt.show_fem == @show_fem... so we need to try to be smart enough to make show_fem not explode
     img = map_img(img, 'resistivity'); % TODO big fat hack to make this work at the expense of an actual function...
     [img, opt] = strip_c2f_background(img, opt, '    ');
     % check we're returning the right size of data
     if isfield(inv_model, 'rec_model')
       img.fwd_model = inv_model.rec_model;
     end
   %  bg = 1;
   %  img.calc_colours.ref_level = bg;
   %  img.calc_colours.clim = bg;
     img.calc_colours.cb_shrink_move = [0.3,0.6,0.02]; % move color bars
     if size(img.elem_data,1) ~= size(img.fwd_model.elems,1)
        warning(sprintf('img.elem_data has %d elements, img.fwd_model.elems has %d elems\n', ...
                        size(img.elem_data,1), ...
                        size(img.fwd_model.elems,1)));
     end
  else % do the "final clean up", same as when we quit
     img = map_img(img, opt.elem_output);
     [img, opt] = strip_c2f_background(img, opt, '    ');
     if isfield(inv_model, 'rec_model')
       img.fwd_model = inv_model.rec_model;
     end
  end
  clf; feval(opt.show_fem, img, 1);
  title(sprintf('iter=%d',k));
  drawnow;
  if isfield(opt,'fig_prefix')
     print('-dpdf',sprintf('%s-fem%d',opt.fig_prefix,k));
     print('-dpng',sprintf('%s-fem%d',opt.fig_prefix,k));
     saveas(gcf,sprintf('%s-fem%d.fig',opt.fig_prefix,k));
  end

% TODO confirm that GN line_search_onm2 is using this residual calculation (preferably, directly)
function [ residual meas elem ] = GN_residual(dv, de, W, hp2RtR)
%   [size(dv); size(W); size(de); size(hp2RtR)]
   % we operate on whatever the iterations operate on (log data, resistance, etc) + perturb(i)*dx
   meas = 0.5*( dv'*W*dv);
   elem = 0.5*(de'*hp2RtR*de);
   residual = meas + elem;

function residual = meas_residual(dv, de, W, hp2RtR)
   residual = norm(dv);

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
%   img = data_mapper(img);
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
%%   img = data_mapper(img,1);
%
%function [img opt] = update_step(org, next, dx, fmin,res, opt)
%   if isfield(opt, 'update_func')
%      [img opt] = feval(opt.update_func,org,next,dx,fmin,res,opt);
%   else
%      img = next;
%   end
%
% function bg = calc_background_resistivity(fmdl, va)
%   % have a look at what we've created
%   % compare data to homgeneous (how good is the model?)
%   % NOTE background conductivity is set by matching amplitude of
%   % homogeneous data against the measurements to get rough matching
%   if(opt.verbose>1)
%     fprintf('est. background resistivity\n');
%   end
%   cache_obj = { fmdl, va };
%   BACKGROUND_R = eidors_obj('get-cache', cache_obj, 'calc_background_resistivity');
%   if isempty(BACKGROUND_R);
%     imgh = mk_image(fmdl, 1); % conductivity = 1 S/m
%     vh = fwd_solve(imgh);
%     % take the best fit of the data
%     BACKGROUND_R = vh.meas \ va; % 32 Ohm.m ... agrees w/ Wilkinson's papers
%     % update cache
%     eidors_obj('set-cache', cache_obj, 'calc_background_resistivity', BACKGROUND_R);
%   else
%     if(opt.verbose > 1)
%       fprintf('  ... cache hit\n');
%     end
%   end
%   if(opt.verbose > 1)
%     fprintf('estimated background resistivity: %0.1f Ohm.m\n', BACKGROUND_R);
%   end

function imdl = deprecate_imdl_opt(imdl,opt)
   if ~isfield(imdl, opt)
      return;
   end
   if ~isstruct(imdl.(opt))
      error(['unexpected inv_model.' opt ' where ' opt ' is not a struct... i do not know what to do']);
   end

   % warn on anything but inv_model.inv_solve.calc_solution_error
   Af = fieldnames(imdl.(opt));
   if ~strcmp(opt, 'inv_solve') || (length(Af(:)) ~= 1) || ~strcmp(Af(:),'calc_solution_error')
      disp(imdl)
      disp(imdl.(opt))
      warning('EIDORS:deprecatedParameters',['INV_SOLVE inv_model.' opt '.* are deprecated in favor of inv_model.inv_solve_abs_core.* as of 30-Apr-2014.']);
   end

   if ~isfield(imdl, 'inv_solve_abs_core')
      imdl.inv_solve_abs_core = imdl.(opt);
   else % we merge
      % merge struct trick from:
      %  http://stackoverflow.com/questions/38645
      for i = fieldnames(imdl.(opt))'
         imdl.inv_solve_abs_core.(i{1})=imdl.(opt).(i{1});
      end
   end
   imdl = rmfield(imdl, opt);

function opt = parse_options(imdl)
   % merge legacy options locations
%   imdl = deprecate_imdl_opt(imdl, 'parameters');
%   imdl = deprecate_imdl_opt(imdl, 'inv_solve');

   % for any general options
   if isfield(imdl, 'inv_solve_abs_core')
      opt = imdl.inv_solve_abs_core;
   else
      opt = struct;
   end

   % verbosity, debug output
   % 0: quiet
   % 1: print iteration count
   % 2: print details as the algorithm progresses
   if ~isfield(opt,'verbose')
      opt.verbose = 4;
      fprintf('  inv_model.inv_solve_abs_core.verbosity not set; defaulting to verbosity=4. See help for details.\n');
   end
   if opt.verbose > 1
      fprintf('  setting default parameters\n');
   end
   % we track how many fwd_solves we do since they are the most expensive part of the iterations
   opt.fwd_solutions = 0;

   if ~isfield(opt, 'show_fem')
      opt.show_fem = @show_fem;
   end

   if ~isfield(opt, 'residual_func') % the objective function
      opt.residual_func = @GN_residual; % r = f(dv, de, W, hp2RtR)
      % NOTE: the meas_residual function exists to maintain
      % compatibility with Nolwenn's code, the GN_residual
      % is a better choice
      %opt.residual_func = @meas_residual; % [r,m,e] = f(dv, de, W, hp2RtR)
   end

   % calculation of update components
   if ~isfield(opt, 'update_func')
      opt.update_func = @GN_update; % dx = f(J, W, hp2RtR, dv, de, opt)
   end
   % figure out if things need to be calculated
   if ~isfield(opt, 'calc_meas_icov') % derivative of the objective function
      opt.calc_meas_icov = 0; % W
   end
   if ~isfield(opt, 'calc_RtR_prior') % derivative of the objective function
      opt.calc_RtR_prior = 0; % RtR
   end
   if ~isfield(opt, 'calc_hyperparameter')
      opt.calc_hyperparameter = 0; % hp2
   end

%   try
      if opt.verbose > 1
         fprintf('    examining function %s(...) for required arguments\n', func2str(opt.update_func));
      end
      % ensure that necessary components are calculated
      % opt.update_func: dx = f(J, W, hp2RtR, dv, de, opt)
%TODO BROKEN      args = function_depends_upon(opt.update_func, 6);
      args = ones(4,1); % TODO BROKEN
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
      disp('    inv_model.inv_solve_abs_core.plot_residuals=0');
   end

   % line search
   if ~isfield(opt,'line_search_func')
      % [alpha, img, dv, opt] = f(img, sx, data0, img0, N, W, hp2RtR, dv, opt);
      opt.line_search_func = @line_search_onm2;
   end
   if ~isfield(opt,'line_search_dv_func')
      opt.line_search_dv_func = @update_dv_core;
      % [dv, opt] = update_dv_core(img, data0, N, opt)
   end
   if ~isfield(opt,'line_search_de_func')
      % we create an anonymous function to skip the first input argument since
      % we always want to calculate de in the line search
      opt.line_search_de_func = @(img, img0, opt) update_de(1, img, img0, opt);
      % de = f(img, img0, opt)
   end
   % an initial guess for the line search step sizes, may be modified by line search
   % TODO this 'sensible default' should be moved to the line_search code since it is not generic to any other line searches
   if ~isfield(opt,'line_search_args') || ...
      ~isfield(opt.line_search_args, 'perturb')
      fmin = 1/4; % arbitrary starting guess
      opt.line_search_args.perturb = [0 fmin/4 fmin/2 fmin fmin*2 fmin*4];
      %opt.line_search_args.perturb = [0 fmin/4 fmin fmin*4];
      %opt.line_search_args.perturb = [0 0.1 0.35 0.7 1.0];
      %opt.line_search_args.perturb = [0 0.1 0.5 0.7 1.0];
      %pt.line_search_args.perturb = [0 0.1 0.7 0.9 1.0];
      %opt.line_search_args.perturb = [0 0.1 0.9 1.0];
   end
   % provide a graphical display of the line search values & fit
   if ~isfield(opt,'line_search_args') || ...
      ~isfield(opt.line_search_args, 'plot')
      if opt.verbose >= 5
         opt.line_search_args.plot = 1;
      else
         opt.line_search_args.plot = 0;
      end
   end
   % pass fig_prefix to the line search as well, unless they are supposed to go somewhere elese
   if isfield(opt,'fig_prefix') && ...
      isfield(opt,'line_search_args') && ...
      ~isfield(opt.line_search_args, 'fig_prefix')
      opt.line_search_args.fig_prefix = opt.fig_prefix;
   end
   % some help on how to turn off the line search plots if we don't want to see them
   if opt.line_search_args.plot ~= 0
      disp('  line search plots (per iteration) are enabled, to disable them set');
      disp('    inv_model.inv_solve_abs_core.line_search_args.plot=0');
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
   if ~isfield(opt, 'c2f_background_fixed')
      opt.c2f_background_fixed = 1; % generally, don't touch the background
   end


   % DATA CONVERSION settings
   % elem type for the initial estimate is based on calc_jacobian_bkgnd which returns an img
   if ~isfield(opt, 'elem_working')
      opt.elem_working = {'conductivity'};
   end
   if ~isfield(opt, 'elem_prior')
      opt.elem_prior = {'conductivity'};
   end
   if ~isfield(opt, 'elem_output')
      opt.elem_output = {'conductivity'};
   end
   if ~isfield(opt, 'meas_input')
      opt.meas_input = 'voltage';
   end
   if ~isfield(opt, 'meas_working')
      opt.meas_working = 'voltage';
   end
   % if the user didn't put these into cell arrays, do
   % so here so there is less error checking later in
   % the code
   for i = {'elem_working', 'elem_prior', 'elem_output'}; %, 'meas_input', 'meas_working'}
     % MATLAB voodoo: deincapsulate a cell containing a
     % string, then use that to access a struct eleemnt
     x = opt.(i{1});
     if ~iscell(x)
        opt.(i{1}) = {x};
     end
   end

   if ~isfield(opt, 'prior_data')
      if isfield(imdl, 'jacobian_bkgnd') && ...
         isfield(imdl.jacobian_bkgnd, 'value') && ...
         length(opt.elem_prior) == 1
         opt.prior_data = {imdl.jacobian_bkgnd.value};
      else
         error('requires inv_model.inv_solve_abs_core.prior_data');
      end
   end

   if ~isfield(opt, 'elem_len')
      if length(opt.elem_working) == 1
         if isfield(imdl.fwd_model, 'coarse2fine')
            c2f = imdl.fwd_model.coarse2fine; % coarse-to-fine mesh mapping
            opt.elem_len = { size(c2f,2) };
         else
            opt.elem_len = { size(imdl.fwd_model.elems,1) };
         end
      else
        error('requires inv_model.inv_solve_abs_core.elem_len');
      end
   end

   % meas_select already handles selecting from the valid measurements
   % we want the same for the elem_data, so we only work on modifying the legal values
   % Note that c2f_background's elements are added to this list if opt.c2f_background_fixed == 1
   if ~isfield(opt, 'elem_fixed') % give a safe default, if none has been provided
      opt.elem_fixed = [];
   elseif iscell(opt.elem_fixed) % if its a cell-array, we convert it to absolute
     % numbers in elem_data
     %  -- requires: opt.elem_len to already be constructed if it was missing
      offset=0;
      ef=[];
      for i=1:length(opt.elem_fixed)
         ef = [ef, opt.elem_fixed{i} + offset];
         offset = offset + opt.elem_len{i};
      end
      opt.elem_fixed = ef;
   end

   % allow a cell array of jacobians
   if ~isfield(opt, 'jacobian')
      opt.jacobian = imdl.fwd_model.jacobian;
   elseif isfield(imdl.fwd_model, 'jacobian')
      imdl.fwd_model
      imdl
      error('inv_model.fwd_model.jacobian and inv_model.inv_solve_abs_core.jacobian should not both exist: it''s ambiguous');
   end
   % defaul hyperparameter is 1
   if ~isfield(opt, 'hyperparameter')
      opt.hyperparameter = {[]};
      for i=1:length(opt.elem_working)
         opt.hyperparameter{i} = 1;
      end
   end
   % if the user didn't put these into cell arrays, do
   % so here so there is less error checking later in
   % the code
   for i = {'elem_len', 'prior_data', 'jacobian', 'hyperparameter'}
     % MATLAB voodoo: deincapsulate a cell containing a
     % string, then use that to access a struct eleemnt
     x = opt.(i{1});
     if ~iscell(x)
        opt.(i{1}) = {x};
     end
   end
   % show what the hyperparameters are configured to when logging
   if opt.verbose > 1
      fprintf('  hyperparameters\n');
      try hp_global = imdl.hyperparameter.value;
          hp_global_str = sprintf(' x %0.4g',hp_global);
      catch hp_global = 1;
          hp_global_str = '';
      end
      for i=1:length(opt.elem_working)
         if isnumeric(opt.hyperparameter{i}) && length(opt.hyperparameter{i}) == 1
            fprintf('    %s: %0.4g\n',opt.elem_working{i}, opt.hyperparameter{i}*hp_global);
         elseif isa(opt.hyperparameter{i}, 'function_handle')
            fprintf('    %s: @%s%s\n',opt.elem_working{i}, func2str(opt.hyperparameter{i}), hp_global_str);
         elseif ischar(opt.hyperparameter{i})
            fprintf('    %s: @%s%s\n',opt.elem_working{i}, opt.hyperparameter{i}, hp_global_str);
         else
            fprintf('    %s: ...\n',opt.elem_working{i});
         end
      end
   end

   % REGULARIZATION RtR
   % for constructing the blockwise RtR matrix
   % can be: explicit matrix, blockwise matrix diagonal, or full blockwise matrix
   % blockwise matrices can be function ptrs or explicit
   if ~isfield(opt, 'RtR_prior')
      if isfield(imdl, 'RtR_prior')
         opt.RtR_prior = {imdl.RtR_prior};
      else
         opt.RtR_prior = {[]}; % null matrix (all zeros)
         warning('missing imdl.inv_solve_abs_core.RtR_prior or imdl.RtR_prior: assuming NO regularization RtR=0');
      end
   end
   % bit of a make work project but if its actually a full numeric matrix we
   % canoncialize it by breaking it up into the blockwise components
   if isnumeric(opt.RtR_prior)
      if size(opt.RtR_prior, 1) ~= size(opt.RtR_prior, 2)
         error('expecting square matrix for imdl.RtR_prior or imdl.inv_solve_abs_core.RtR_prior');
      end
      if length(opt.RtR_prior) == 1
         opt.RtR_prior = {opt.RtR_prior}; % encapsulate directly into a cell array
      else
         RtR = opt.RtR_prior;
         opt.RtR_prior = {[]};
         esi = 0; eei = 0;
         for i=1:length(opt.elem_len)
            esi = eei +1;
            eei = eei +opt.elem_len{i};
            esj = 0; eej = 0;
            for j=1:length(opt.elem_len)
               esj = eej +1;
               eej = eej +opt.elem_len{j};
               opt.RtR_prior(i,j) = RtR(esi:eei, esj:eej);
            end
         end
      end
   elseif ~iscell(opt.RtR_prior) % not a cell array? encapsulate it
      opt.RtR_prior = {opt.RtR_prior};
   end
   % if not square then expand the block matrix
   % single row/column: this is our diagonal --> expand to full blockwise matrix
   if any(size(opt.RtR_prior) ~= ([1 1]*length(opt.elem_len)))
      if (size(opt.RtR_prior, 1) ~= 1) && ...
         (size(opt.RtR_prior, 2) ~= 1)
         error('odd imdl.RtR_prior or imdl.inv_solve_abs_core.RtR_prior, cannot figure out how to expand it blockwise');
      end
      if (size(opt.RtR_prior, 1) ~= length(opt.elem_len)) && ...
         (size(opt.RtR_prior, 2) ~= length(opt.elem_len))
         error('odd imdl.RtR_prior or imdl.inv_solve_abs_core.RtR_prior, not enough blockwise components vs. elem_working types');
      end
      RtR_diag = opt.RtR_prior;
      opt.RtR_prior = {[]}; % delete and start again
      for i=1:length(opt.elem_len)
         opt.RtR_prior(i,i) = RtR_diag(i);
      end
   end

   % now sort out the hyperparameter for the "R^T R" (RtR) matrix
   hp=opt.hyperparameter;
   if size(hp,2) == 1 % one column
      hp = hp'; % ... now one row
   end
   if iscell(hp)
      % if it's a cell array that matches size of the RtR, then we're done
      if all(size(hp) == size(opt.RtR_prior))
         opt.hyperparameter = hp;
      % if the columns matches, then we can expand on the diangonal, everything else gets '1'
      elseif length(hp) == length(opt.RtR_prior)
         opt.hyperparameter = opt.RtR_prior;
         [opt.hyperparameter{:}] = deal(1); % hp = 1 everywhere
         opt.hyperparameter(logical(eye(size(opt.RtR_prior)))) = hp; % assign to diagonal
      else
         error('hmm, don''t understand this opt.hyperparameter cellarray');
      end
   % if it's a single hyperparameter, that's the value everywhere
   elseif isnumeric(hp)
      opt.hyperparameter = opt.RtR_prior;
      [opt.hyperparameter{:}] = deal({hp});
   else
      error('don''t understand this opt.hyperparameter');
   end

   % JACOBIAN CHAIN RULE conductivity -> whatever
   % where x = conductivity at this iteration
   %       S = a scaling matrix, generally a diagonal matrix of size matching Jacobian columns
   % Jn = J * S;
   % if not provided, determine based on 'elem_working' type
   if ~isfield(opt, 'calc_jacobian_scaling_func')
      pinv = strfind(opt.elem_working, 'resistivity');
      plog = strfind(opt.elem_working, 'log_');
      plog10 = strfind(opt.elem_working, 'log10_');
      for i = 1:length(opt.elem_working)
        prefix = '';
        if plog{i}
           prefix = 'log';
        elseif plog10{i}
           prefix = 'log10';
        else
           prefix = '';
        end
        if pinv{i}
           prefix = [prefix '_inv'];
        end
        switch(prefix)
          case ''
             opt.calc_jacobian_scaling_func{i} = @ret1_func;  % S = f(x)
          case 'log'
             opt.calc_jacobian_scaling_func{i} = @dx_dlogx;   % S = f(x)
          case 'log10'
             opt.calc_jacobian_scaling_func{i} = @dx_dlog10x; % S = f(x)
          case '_inv'
             opt.calc_jacobian_scaling_func{i} = @dx_dy;      % S = f(x)
          case 'log_inv'
             opt.calc_jacobian_scaling_func{i} = @dx_dlogy;   % S = f(x)
          case 'log10_inv'
             opt.calc_jacobian_scaling_func{i} = @dx_dlog10y; % S = f(x)
          otherwise
             error('oops');
       end
     end
   end

   if ~isfield(opt, 'update_img_func')
      opt.update_img_func = @null_func; % img = f(img, opt)
   end

   if ~isfield(opt, 'return_working_variables')
      opt.return_working_variables = 0;
   end

function check_matrix_sizes(J, W, hp2RtR, dv, de, opt)
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
   if any(size(J) ~= [nv ne])
      error('J size (%d rows, %d cols) is incorrect (%d rows, %d cols)', size(J), nv, ne);
   end
   if opt.calc_RtR_prior && ...
      any(size(hp2RtR) ~= [ne ne])
      error('hp2RtR size (%d rows, %d cols) is incorrect (%d rows, %d cols)', size(hp2RtR), ne, ne);
   end

function dx = update_dx(J, W, hp2RtR, dv, de, opt)
   if(opt.verbose > 1)
      fprintf( '    calc step size dx\n');
   end

   % don't penalize for fixed elements
   de(opt.elem_fixed) = 0;

   % TODO move this outside the inner loop of the iterations, it only needs to be done once
   check_matrix_sizes(J, W, hp2RtR, dv, de, opt)

   % zero out the appropriate things so that we can get a dx=0 for the elem_fixed
   J(:,opt.elem_fixed) = 0;
   de(opt.elem_fixed,:) = 0;
   hp2RtR(opt.elem_fixed,:) = 0;
   V=opt.elem_fixed;
   N=size(hp2RtR,1)+1;
   hp2RtR(N*(V-1)+1) = 1; % set diagonals to 1 to avoid divide by zero
   % do the update step direction calculation
   dx = feval(opt.update_func, J, W, hp2RtR, dv, de, opt);

   % check that our elem_fixed stayed fixed
   if any(dx(opt.elem_fixed) ~= 0)
      error('elem_fixed did''t give dx=0 at update_dx')
   end

   if(opt.verbose > 1)
      fprintf('      ||dx||=%0.3g\n', norm(dx));
      es = 0; ee = 0;
      for i=1:length(opt.elem_working)
          es = ee +1; ee = ee + opt.elem_len{i};
          nd = norm(dx(es:ee));
          fprintf( '      ||dx_%d||=%0.3g (%s)\n',i, nd, opt.elem_working{i});
      end
   end

function dx = GN_update(J, W, hp2RtR, dv, de, opt)
   try
      % the actual update
      dx = (J'*W*J + hp2RtR)\(J'*dv + hp2RtR*de); % LU decomp
   catch ME % boom
      tol = 1e-6; % default 1e-6
      maxit = []; % default [] --> min(n,20)
      M = []; % default [] --> no preconditioner
      x0 = []; % default [] --> zeros(n,1)

      % try Preconditioned Conjugate Gradient: A x = b, solve for x
      % avoids J'*J for n x m matrix with large number of m cols --> J'*J becomes an m x m dense matrix
      LHS = @(x) J'*(W*(J*x)) + hp2RtR*x;
      RHS = J'*dv + hp2RtR*de;

      tol=100*eps*size(J,2)^2; % rough estimate based on multiply-accumulates
%      maxit = 10;
      [dx, flag, relres, iter, resvec] = pcg(LHS, RHS, tol, maxit, M, x0);
      % TODO if verbose...
      switch flag
         case 0
            if opt.verbose > 1
               fprintf('      PCG: relres=%g < tol=%g @ iter#%d\n',relres,tol,iter);
            end
         case 1
            if opt.verbose > 1
               fprintf('      PCG: relres=%g > tol=%g @ iter#%d (max#%d) [max iter]\n',relres,tol,iter,maxit);
            end
         case 2
            error('error: PCG ill-conditioned preconditioner M');
         case 3
            if opt.verbose > 1
               fprintf('      PCG: relres=%g > tol=%g @ iter#%d (max#%d) [stagnated]\n',relres,tol,iter,maxit);
            end
         case 4
            error('error: PCG a scalar quantity became too large or small to continue');
         otherwise
            error(sprintf('error: PCG unrecognized flag=%d',flag));
      end
% NOTE PCG is still a work in progress and generally problem specific
%      % plot convergence for pcg()
%      clf;
%         xlabel('iteration #');
%         ylabel('relative residual');
%         xx = 0:length(resvec)-1;
%         semilogy(xx,resvec/norm(RHS),'b.');
%         hold on;
%         legend('no preconditioner');
%hold on
%% sparsity strategies: www.cerfacs.fr/algor/reports/Dissertations/TH_PA_02_48.pdf
%sJ = J; sJ(J/max(J(:)) > 0.005) = 0; sJ=sparse(sJ); % sparsify J
%M = ichol(sJ'*W*sJ + hp2RtR + speye(length(J)));
%      [dx, flag, relres, iter, resvec] = pcg(LHS, RHS, tol, maxit, M);
%         semilogy(xx,resvec/norm(RHS),'r.');
%         legend('no P', 'IC(sp(J'')*W*sp(J) + hp2RtR)');
%
%
%clf; Jj=J(:,1:800); imagesc(Jj'*Jj);

      % compare with gmres, bicgstab, lsqr
      % try preconditioners ilu, ichol (incomplete LU or Cholesky decomp.)

      %rethrow(ME); % we assume this is an 'excessive memory requested' failure
   end

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
function out = null_func(in, varargin);
   out = in;

% this function always returns one
function [out, x, y, z] = ret1_func(varargin);
   out = 1;
   x = [];
   y = [];
   z = [];

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
      if opt.c2f_background_fixed
         % TODO assumes conductivity/resistivity is the *first* parameterization
         opt.elem_fixed(end+1) = nc+1;
      end
    end
    % TODO assumes conductivity/resistivity is the *first* parameterization
    opt.elem_len(1) = {size(c2f,2)}; % elem_len +1

function [img, opt] = strip_c2f_background(img, opt, indent)
    if nargin < 3
       indent = '';
    end
    % nothing to do?
    if opt.c2f_background <= 0
      return;
    end

    % if there are multiple 'params' (parametrizations), we assume its the first
    % TODO -- this isn't a great assumption but it'll work for now,
    %         we should add a better (more general) mechanism
    in = img.current_params;
    out = opt.elem_output;
    if iscell(in)
       in = in{1};
    end
    if iscell(out)
       out = out{1};
    end

    % go about cleaning up the background
    e = opt.c2f_background;
    % take backgtround elements and convert to output
    % 'params' (resistivity, etc)
    bg = map_data(img.elem_data(e), in, out);
    img.elem_data_background = bg;
    % remove elements from elem_data & c2f
    img.elem_data(e) = [];
    img.fwd_model.coarse2fine(:,e) = [];
    % remove our element from the lists
    opt.c2f_background = 0;
    ri = find(opt.elem_fixed == e);
    opt.elem_fixed(ri) = [];
    if isfield(img, 'params_sel')
       for i = 1:length(img.params_sel)
          t = img.params_sel{i};
          ti = find(t == e);
          t(ti) = []; % rm 'e' from the list of params_sel
          ti = find(t > e);
          t(ti) = t(ti)-1; % down-count element indices greater than our deleted one
          img.params_sel{i} = t;
       end
    end

    % show what we got for a background value
    if(opt.verbose > 1)
       bg = img.elem_data_background;
       bg = map_data(bg, in, 'resistivity');
       fprintf('%s  background conductivity: %0.1f Ohm.m\n', indent, bg);
    end

function b = has_params(s)
b = false;
if isstruct(s)
   b = any(ismember(fieldnames(s),supported_params));
end

% wrapper function for to_base_types
function out = map_img_base_types(img)
  out = to_base_types(img.current_params);

% convert from know types to their base types
% A helper function for getting to a basic paramterization
% prior to any required scaling, etc.
function type = to_base_types(type)
  if ~iscell(type)
     type = {type};
  end
  for i = 1:length(type);
     type(i) = {strrep(type{i}, 'log_', '')};
     type(i) = {strrep(type{i}, 'log10_', '')};
     type(i) = {strrep(type{i}, 'resistivity', 'conductivity')};
     type(i) = {strrep(type{i}, 'apparent_resistivity', 'voltage')};
  end

function img = map_img(img, out);
   err_if_inf_or_nan(img.elem_data, 'img-pre');
   try in = img.current_params;
   catch in = {'conductivity'};
   end
   % make cell array of strings
   if ischar(in)
      in = {in};
      img.current_params = in;
   end
   if ischar(out)
      out = {out};
   end

   % if we have mixed data, check that we have a selector to differentiate between them
   if ~isfield(img, 'params_sel')
      if length(in(:)) == 1
         img.params_sel = {1:size(img.elem_data,1)};
      else
         error('found multiple parametrizations (params) but no params_sel cell array in img');
      end
   end

   % create data?! we don't know how
   if length(out(:)) > length(in(:))
      error('missing data (more out types than in types)');
   elseif length(out(:)) < length(in(:))
      % delete data: we can do that
      % TODO we could genearlize this into a reorganizing tool BUT we're just
      % interested in something that works, so if we have more than one out(:),
      % we don't know what to do currently and error out
      if length(out(:)) ~= 1
         error('map_img can convert ALL parameters or select a SINGLE output type from multiple input types');
      end
      inm  = to_base_types(in);
      outm = to_base_types(out);
      del = sort(find(~strcmp(outm(1), inm(:))), 'descend'); % we do this loop backwards in the hopes of avoiding shuffling data that is about to be deleted
      if length(del)+1 ~= length(in)
         error('Confused about what to remove from the img. You''ll need to sort the parametrizations out yourself when removing data.');
      end
      for i = del(:)' % delete each of the extra indices
         ndel = length(img.params_sel{i}); % number of deleted elements
         for j = i+1:length(img.params_sel)
            img.params_sel{j} = img.params_sel{j} - ndel;
         end
         img.elem_data(img.params_sel{i}) = []; % rm elem_data
         img.params_sel(i) = []; % rm params_sel
         img.current_params(i) = []; % rm current_params
      end
      in = img.current_params;
   end

   % the sizes now match, we can do the mapping
   for i = 1:length(out(:))
      % map the data
      x = img.elem_data(img.params_sel{i});
      x = map_data(x, in{i}, out{i});
      img.elem_data(img.params_sel{i}) = x;
      img.current_params{i} = out{i};
   end
   err_if_inf_or_nan(img.elem_data, 'img-post');

   % clean up params_sel/current_params if we only have one parametrization
   if length(img.current_params(:)) == 1
      img.current_params = img.current_params{1};
      img = rmfield(img, 'params_sel'); % unnecessary since we know its all elem_data
   end

function x = map_data(x, in, out)
   % check that in and out are single strings, not lists of strings
   if ~ischar(in)
      if iscell(in) && (length(in(:)) == 1)
         in = in{1};
      else
         error('expecting single string for map_data() "in" type');
      end
   end
   if ~ischar(out)
      if iscell(out) && (length(out(:)) == 1)
         out = out{1};
      else
         error('expecting single string for map_data() "out" type');
      end
   end

   % quit early if there is nothing to do
   if strcmp(in, out) % in == out
      return; % do nothing
   end

   % resistivity to conductivity conversion
   % we can't get here if in == out
   % we've already checked for log convserions on input or output
   if any(strcmp(in,  {'resistivity', 'conductivity'})) && ...
      any(strcmp(out, {'resistivity', 'conductivity'}))
      x = 1./x; % conductivity <-> resistivity
   % log conversion
   elseif any(strcmp({in(1:3), out(1:3)}, 'log'))
      % log_10 x -> x
      if strcmp(in(1:6), 'log10_')
         if any(x >= log10(realmax)-eps) warning('loss of precision -> inf'); end
         x = map_data(10.^x, in(7:end), out);
      % ln x -> x
      elseif strcmp(in(1:4), 'log_')
         if any(x >= log(realmax)-eps) warning('loss of precision -> inf'); end
         x = map_data(exp(x), in(5:end), out);
      % x -> log_10 x
      elseif strcmp(out(1:6), 'log10_')
         if any(x <= 0 + eps) warning('loss of precision -> -inf'); end
         x = log10(map_data(x, in, out(7:end)));
      % x -> ln x
      elseif strcmp(out(1:4), 'log_')
         if any(x <= 0 + eps) warning('loss of precision -> -inf'); end
         x = log(map_data(x, in, out(5:end)));
      else
         error(sprintf('unknown conversion (log conversion?) %s - > %s', in, out));
      end
   else
      error('unknown conversion %s -> %s', in, out);
   end
   x(x == +inf) = +realmax;
   x(x == -inf) = -realmax;
   err_if_inf_or_nan(x, 'map_data-post');

function b = map_meas(b, N, in, out)
   err_if_inf_or_nan(b, 'map_meas-pre');
   if strcmp(in, out) % in == out
      return; % do nothing
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
   % log conversion
   elseif any(strcmp({in(1:3), out(1:3)}, 'log'))
      % log_10 b -> b
      if strcmp(in(1:6), 'log10_')
         if any(b > log10(realmax)-eps) warning('loss of precision -> inf'); end
         b = map_meas(10.^b, N, in(7:end), out);
      % ln b -> b
      elseif strcmp(in(1:4), 'log_')
         if any(b > log(realmax)-eps) warning('loss of precision -> inf'); end
         b = map_meas(exp(b), N, in(5:end), out);
      % b -> log_10 b
      elseif strcmp(out(1:6), 'log10_')
         if any(b <= 0+eps) warning('loss of precision -> -inf'); end
         b = log10(map_meas(b, N, in, out(7:end)));
      % b -> ln b
      elseif strcmp(out(1:4), 'log_')
         if any(b <= 0+eps) warning('loss of precision -> -inf'); end
         b = log(map_meas(b, N, in, out(5:end)));
      else
         error(sprintf('unknown conversion (log conversion?) %s - > %s', in, out));
      end
   else
      error('unknown conversion %s -> %s', in, out);
   end
   err_if_inf_or_nan(b, 'map_meas-post');

function x=range(y)
x=max(y)-min(y);

function do_unit_test(solver)
   if nargin == 0
     solver = 'inv_solve_abs_core';
   end
   do_unit_test_rec_mv(solver);
   do_unit_test_sub;
   do_unit_test_rec1(solver);
%pass = pass & do_unit_test_rec2(solver);
% TODO the ..._rec2 unit test is very, very slow... what can we do to speed it up... looks like the perturbations get kinda borked when using the line_search_onm2

% test sub-functions
% map_meas, map_data
% jacobian scalings
function do_unit_test_sub
d = 1;
while d ~= 1 & d ~= 0
  d = rand(1);
end
disp('TEST: map_data()');
elem_types = {'conductivity', 'log_conductivity', 'log10_conductivity', ...
              'resistivity',  'log_resistivity',  'log10_resistivity'};
expected = [d         log(d)         log10(d)      1./d      log(1./d)      log10(1./d); ...
            exp(d)    d              log10(exp(d)) 1./exp(d) log(1./exp(d)) log10(1./exp(d)); ...
            10.^d     log(10.^d )    d             1./10.^d  log(1./10.^d ) log10(1./10.^d ); ...
            1./d      log(1./d  )    log10(1./d)   d         log(d)         log10(d); ...
            1./exp(d) log(1./exp(d)) log10(1./exp(d)) exp(d) d              log10(exp(d)); ...
            1./10.^d  log(1./10.^d)  log10(1./10.^d)  10.^d  log(10.^d)     d ];
for i = 1:length(elem_types)
  for j = 1:length(elem_types)
    test_map_data(d, elem_types{i}, elem_types{j}, expected(i,j));
  end
end

disp('TEST: map_meas()');
N = 1/15;
Ninv = 1/N;
% function b = map_meas(b, N, in, out)
elem_types = {'voltage', 'log_voltage', 'log10_voltage', ...
              'apparent_resistivity',  'log_apparent_resistivity',  'log10_apparent_resistivity'};
expected = [d         log(d)         log10(d)      N*d      log(N*d)      log10(N*d); ...
            exp(d)    d              log10(exp(d)) N*exp(d) log(N*exp(d)) log10(N*exp(d)); ...
            10.^d     log(10.^d )    d             N*10.^d  log(N*10.^d ) log10(N*10.^d ); ...
            Ninv*d      log(Ninv*d  )    log10(Ninv*d)   d         log(d)         log10(d); ...
            Ninv*exp(d) log(Ninv*exp(d)) log10(Ninv*exp(d)) exp(d) d              log10(exp(d)); ...
            Ninv*10.^d  log(Ninv*10.^d)  log10(Ninv*10.^d)  10.^d  log(10.^d)     d ];
for i = 1:length(elem_types)
  for j = 1:length(elem_types)
    test_map_meas(d, N, elem_types{i}, elem_types{j}, expected(i,j));
  end
end

disp('TEST: Jacobian scaling');
d = [d d]';
unit_test_cmp( ...
   sprintf('Jacobian scaling (%s)', 'conductivity'), ...
   ret1_func(d), 1);

unit_test_cmp( ...
   sprintf('Jacobian scaling (%s)', 'log_conductivity'), ...
   dx_dlogx(d), diag(d));

unit_test_cmp( ...
   sprintf('Jacobian scaling (%s)', 'log10_conductivity'), ...
   dx_dlog10x(d), diag(d)*log(10));

unit_test_cmp( ...
   sprintf('Jacobian scaling (%s)', 'resistivity'), ...
   dx_dy(d), diag(-d.^2));

unit_test_cmp( ...
   sprintf('Jacobian scaling (%s)', 'log_resistivity'), ...
   dx_dlogy(d), diag(-d));

unit_test_cmp( ...
   sprintf('Jacobian scaling (%s)', 'log10_resistivity'), ...
   dx_dlog10y(d), diag(-d)/log(10));


function test_map_data(data, in, out, expected)
%fprintf('TEST: map_data(%s -> %s)\n', in, out);
   calc_val = map_data(data, in, out);
   str = sprintf('map_data(%s -> %s)', in, out);
   unit_test_cmp(str, calc_val, expected)

function test_map_meas(data, N, in, out, expected)
%fprintf('TEST: map_meas(%s -> %s)\n', in, out);
   calc_val = map_meas(data, N, in, out);
   str = sprintf('map_data(%s -> %s)', in, out);
   unit_test_cmp(str, calc_val, expected)


% a couple easy reconstructions
% check c2f, apparent_resistivity, log_conductivity, verbosity don't error out
function do_unit_test_rec1(solver)
% -------------
% ADAPTED FROM
% Create simulation data $Id$
%  http://eidors3d.sourceforge.net/tutorial/adv_image_reconst/basic_iterative.shtml
% 3D Model
imdl= mk_common_model('c2t4',16); % 576 elements
imdl.solve = solver;
imdl.reconst_type = 'absolute';
imdl.inv_solve_abs_core.prior_data = 1;
imdl.inv_solve_abs_core.elem_prior = 'conductivity';
imdl.inv_solve_abs_core.elem_working = 'log_conductivity';
imdl.inv_solve_abs_core.meas_working = 'apparent_resistivity';
imdl.inv_solve_abs_core.calc_solution_error = 0;
imdl.inv_solve_abs_core.verbose = 0;
%show_fem(imdl.fwd_model);
imgsrc= mk_image( imdl.fwd_model, 1);
% set homogeneous conductivity and simulate
vh=fwd_solve(imgsrc);
% set inhomogeneous conductivity and simulate
ctrs= interp_mesh(imdl.fwd_model);
x= ctrs(:,1); y= ctrs(:,2);
r1=sqrt((x+5).^2 + (y+5).^2); r2 = sqrt((x-85).^2 + (y-65).^2);
imgsrc.elem_data(r1<50)= 0.05;
imgsrc.elem_data(r2<30)= 100;
imgp = map_img(imgsrc, 'log10_conductivity');
hh=clf; subplot(221); show_fem(imgp,1); axis tight; title('synthetic data, logC');
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
imdl.inv_solve_abs_core.verbose = 0;
%imdl.inv_solve_abs_core.meas_working = 'apparent_resistivity';
img2= inv_solve(imdl, vi);
figure(hh); subplot(223); show_fem(img2,1); axis tight; title('#2 verbosity=0');
max_err = max(abs((img1.elem_data - img2.elem_data)./(img1.elem_data)));

unit_test_cmp('img1 == img2', max_err >0.15, 0);
if max_err > 0.15
  fprintf('TEST:  img1 != img2 --> FAIL %g %%\n', max_err);
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
%imdl.inv_solve_abs_core.verbose = 10;
img3= inv_solve(imdl, vi);
%figure(hh); subplot(224); show_fem(cmdl,1); axis tight; title('#3 c2f');
figure(hh); subplot(223); show_fem(img3,1); axis tight; title('#3 c2f');
% check
e1 = c2f \ img1.elem_data; % noisy and unstable... but its a crude check
e3 = img3.elem_data;
err = abs((e1 - e3) ./ e1);
err(abs(e1) < 20) = 0;
err_thres = 0.55;

unit_test_cmp('img1 == img3', any(err > err_thres), 0);
if any(err > err_thres) % maximum 15% error
  ni = find(err > err_thres);
  fprintf('TEST:  img1 != img3 --> FAIL max(err) = %0.2e on %d elements (thres=%0.2e)\n', ...
          max(err(ni)), length(ni), err_thres);
end

%imdl.inv_solve_abs_core.verbose = 1000;
imdl.inv_solve_abs_core.elem_output = 'log10_resistivity'; % resistivity output works
img4= inv_solve(imdl, vi);
figure(hh); subplot(224); show_fem(img4,1); axis tight; title('#4 c2f + log10 resistivity out');
% check
e4 = 1./(10.^img4.elem_data);
err = abs((e1 - e4) ./ e1);
err(abs(e1) < 20) = 0;
err_thres = 0.40;

unit_test_cmp('img1 == img4', any(err > err_thres), 0);
if any(err > err_thres) % maximum 15% error
  ni = find(err > err_thres);
  fprintf('TEST:  img1 != img4 --> FAIL max(err) = %0.2e on %d elements (thres=%0.2e)\n', ...
          max(err(ni)), length(ni), err_thres);
end

% a couple easy reconstructions with movement or similar
function do_unit_test_rec_mv(solver)
disp('TEST: conductivity and movement --> baseline conductivity only');
% -------------
% ADAPTED FROM
% Create simulation data $Id$
%  http://eidors3d.sourceforge.net/tutorial/adv_image_reconst/basic_iterative.shtml
% 3D Model
imdl= mk_common_model('c2t4',16); % 576 elements
ne = length(imdl.fwd_model.electrode);
nt = length(imdl.fwd_model.elems);
imdl.solve = solver;
imdl.reconst_type = 'absolute';
% specify the units to work in
imdl.inv_solve_abs_core.meas_input   = 'voltage';
imdl.inv_solve_abs_core.meas_working = 'apparent_resistivity';
imdl.inv_solve_abs_core.elem_prior   = {   'conductivity'   };
imdl.inv_solve_abs_core.prior_data   = {        1           };
imdl.inv_solve_abs_core.elem_working = {'log_conductivity'};
imdl.inv_solve_abs_core.elem_output  = {'log10_conductivity'};
imdl.inv_solve_abs_core.calc_solution_error = 0;
imdl.inv_solve_abs_core.verbose = 0;
imdl.hyperparameter.value = 0.01;

% set homogeneous conductivity and simulate
imgsrc= mk_image( imdl.fwd_model, 1);
vh=fwd_solve(imgsrc);
% set inhomogeneous conductivity
ctrs= interp_mesh(imdl.fwd_model);
x= ctrs(:,1); y= ctrs(:,2);
r1=sqrt((x+5).^2 + (y+5).^2); r2 = sqrt((x-45).^2 + (y-35).^2);
imgsrc.elem_data(r1<50)= 0.05;
imgsrc.elem_data(r2<30)= 100;

% inhomogeneous data
vi=fwd_solve( imgsrc );
% add noise
%Add 30dB SNR noise to data
noise_level= std(vi.meas - vh.meas)/10^(30/20);
vi.meas = vi.meas + noise_level*randn(size(vi.meas));

% show model
hh=clf; subplot(221); imgp = map_img(imgsrc, 'log10_conductivity'); show_fem(imgp,1); axis tight; title('synth baseline, logC');

% Reconstruct Images
img0= inv_solve(imdl, vi);
figure(hh); subplot(222);
 img0 = map_img(img0, 'log10_conductivity');
 show_fem(img0, 1); axis tight;

disp('TEST: conductivity + movement');
imdl.fwd_model = rmfield(imdl.fwd_model, 'jacobian');
% specify the units to work in
imdl.inv_solve_abs_core.elem_prior   = {   'conductivity'   , 'movement'};
imdl.inv_solve_abs_core.prior_data   = {        1           ,     0     };
imdl.inv_solve_abs_core.RtR_prior    = {     @eidors_default, @prior_movement_only};
imdl.inv_solve_abs_core.elem_len     = {       nt           ,   ne*2    };
imdl.inv_solve_abs_core.elem_working = {  'log_conductivity', 'movement'};
imdl.inv_solve_abs_core.elem_output  = {'log10_conductivity', 'movement'};
imdl.inv_solve_abs_core.jacobian     = { @jacobian_adjoint  , @jacobian_movement_only};
imdl.inv_solve_abs_core.hyperparameter = {   [1 1.1 0.9]    ,  sqrt(2e-3)     }; % multiplied by imdl.hyperparameter.value
imdl.inv_solve_abs_core.verbose = 2;

% electrode positions before
nn = [imgsrc.fwd_model.electrode(:).nodes];
elec_orig = imgsrc.fwd_model.nodes(nn,:);
% set 2% radial movement
nn = imgsrc.fwd_model.nodes;
imgsrc.fwd_model.nodes = nn * [1-0.02 0; 0 1+0.02]; % 1% compress X, 1% expansion Y, NOT conformal
% electrode positions after
nn = [imgsrc.fwd_model.electrode(:).nodes];
elec_mv = imgsrc.fwd_model.nodes(nn,:);

% inhomogeneous data
vi=fwd_solve( imgsrc );
% add noise
%Add 30dB SNR noise to data
noise_level= std(vi.meas - vh.meas)/10^(30/20);
%vi.meas = vi.meas + noise_level*randn(size(vi.meas));

% show model
nn = [imgsrc.fwd_model.electrode(1:4).nodes];
figure(hh); subplot(223); imgp = map_img(imgsrc, 'log10_conductivity'); show_fem_move(imgp,elec_mv-elec_orig,10,1); axis tight; title('synth mvmt, logC');

% Reconstruct Images
img1= inv_solve(imdl, vi);
figure(hh); subplot(224);
 imgm = map_img(img1, 'movement');
 img1 = map_img(img1, 'log10_conductivity');
 show_fem_move(img1,reshape(imgm.elem_data,16,2), 10, 1); axis tight;

% TEST for mismatch on coductivity image
err = abs((img0.elem_data - img1.elem_data) ./ img0.elem_data);
err(abs(img0.elem_data)/max(abs(img0.elem_data)) < 0.50) = 0;
err_thres = 0.40;

unit_test_cmp('img0 == img1 + mvmt', any(err > err_thres), 0);

if any(err > err_thres) % maximum 15% error
  ni = find(err > err_thres);
  fprintf('TEST:  img0 != img1 + mvmt --> FAIL max(err) = %0.2e on %d elements (thres=%0.2e)\n', ...
          max(err(ni)), length(ni), err_thres);
end

% helper function: calculate jacobian movement by itself
function Jm = jacobian_movement_only (fwd_model, img);
  pp = fwd_model_parameters(img.fwd_model);
  szJm = pp.n_elec * pp.n_dims; % number of electrodes * dimensions
  img = map_img(img, 'conductivity'); % expect conductivity only
  Jcm = jacobian_movement(fwd_model, img);
  Jm = Jcm(:,(end-szJm+1):end);
%% this plot shows we are grabing the right section of the Jacobian
%  figure();
%  subplot(311); imagesc(Jcm); axis ij equal tight; xlabel(sprintf('||Jcm||=%g',norm(Jcm))); colorbar;
%  Jc = jacobian_adjoint(fwd_model, img);
%  subplot(312); imagesc([Jc Jm]); axis ij equal tight; xlabel(sprintf('||[Jc Jm]||=%g',norm([Jc Jm]))); colorbar;
%  dd = abs([Jc Jm]-Jcm); % difference
%  subplot(313); imagesc(dd); axis ij equal tight; xlabel(sprintf('|| |[Jc Jm]-Jcm| ||=%g',norm(dd))); colorbar;

function RtR = prior_movement_only(imdl);
  imdl.image_prior.parameters(1) = 1; % weighting of movement vs. conductivity ... but we're dropping conductivity here
  pp = fwd_model_parameters(imdl.fwd_model);
  szPm = pp.n_elec * pp.n_dims; % number of electrodes * dimensions
  RtR = prior_movement(imdl);
  RtR = RtR((end-szPm+1):end,(end-szPm+1):end);

function do_unit_test_rec2(solver)
disp('TEST: reconstruct a discontinuity');
shape_str = ['solid top    = plane(0,0,0;0,1,0);\n' ...
             'solid mainobj= top and orthobrick(-100,-200,-100;410,10,100) -maxh=20.0;\n'];
e0 = linspace(0,310,64)';
elec_pos = [e0,0*e0,0*e0,1+0*e0,0*e0,0*e0];
elec_shape= [0.1,0.1,1];
elec_obj = 'top';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);
%fmdl.nodes = fmdl.nodes(:,[1,3,2]);
% spacing= [1 1 1 2 3 3 4 4 5 6 6 7 8 8 9 10 10 11 12 12 13 14 14 15 16 17];
% multiples= [1 2 3 2 1 5/3 1 2  1 1 7/6 1 1 10/8 1 1 12/10 1 1 13/12 1 1 15/14 1 1 1];
% fmdl.stimulation= stim_pattern_geophys( 64, 'Schlumberger', {'spacings', spacing,'multiples',multiples});

fmdl.stimulation= stim_pattern_geophys( 64, 'Wenner', {'spacings', 1:32} );

cmdl= mk_grid_model([], 2.5+[-30,5,20,30:10:290,300,315,340], ...
                            -[0:5:10 17 30 50 75 100]);
% having a c2f on the coarse model f#$%s up the c2f calculator
cmdl= rmfield(cmdl, 'coarse2fine');
% cmdl = mk_grid_model([], 2.5+[-50,-20,0:10:310,330,360], ...
%                              -[0:2.5:10, 15:5:25,30:10:80,100,120]);
[c2f, b_c2f] = mk_coarse_fine_mapping( fmdl, cmdl);
% c2f maps cmdl elements to fmdl elements, where there are no cmdl elements
% b_c2f is the cmdl background element that is mapped to the fmdl elements
% Note: adding a background element, this is now done inside the inv_solve_abs_GN solver if required
%S= sum(c2f,2);
% find fractional c2f elements
%b= find(S<0.9999);
% find almost complete c2f elements
%a= find(S>=0.9999 & S<1);
% fix potential rounding problems by normalizing c2f to sum to 1
%c2f(a,:)= c2f(a,:)./repmat(S(a),1,size(c2f,2));
% remove the entire element's mapping and assign it the background conductivity
%c2f(b,:)= 0; c2f(b,end+1)= 1;
fmdl.coarse2fine= c2f;

% generate sythetic data
img = mk_image(fmdl,1);
fm_pts = interp_mesh(fmdl);
x_bary= fm_pts(:,1); z_bary= fm_pts(:,2);
z_params= (min(fmdl.nodes(:,2)):max(fmdl.nodes(:,2)))';
a = 0.36;
b = 130;
x_params= a*z_params+b;
xlim=interp1(z_params,x_params,z_bary);
img.elem_data(x_bary>xlim)= 0.01;
clf; show_fem(img); title('model');

% img2= mk_image(fmdl,img.elem_data);
% clf; show_fem(img2);

% img = mk_image(fmdl,0+ mk_c2f_circ_mapping(fmdl,[100;-30;0;50])*100);
% img.elem_data(img.elem_data==0)= 0.1;
dd  = fwd_solve(img);
% TODO add some noise!!!

imdl= eidors_obj('inv_model','test');
imdl.fwd_model= fmdl;
imdl.rec_model= cmdl;
imdl.fwd_model.normalize_measurements = 0;
imdl.rec_model.normalize_measurements = 0;
imdl.RtR_prior = @prior_laplace;
%imdl.RtR_prior = @prior_tikhonov;
imdl.solve = solver;
imdl.reconst_type = 'absolute';
imdl.hyperparameter.value = 1e2; % was 0.1
imdl.jacobian_bkgnd.value = 1;

imdl.inv_solve_abs_core.elem_working = 'log_conductivity';
imdl.inv_solve_abs_core.meas_working = 'apparent_resistivity';
imdl.inv_solve_abs_core.dtol_iter = 4; % default 1 -> start checking on the first iter
imdl.inv_solve_abs_core.max_iterations = 20; % default 10

% the conversion to apparaent resistivity is now handled inside the solver
%%img1= mk_image(fmdl,1);
%%vh1= fwd_solve(img1);
%%normalisation= 1./vh1.meas;
%%I= speye(length(normalisation));
%%I(1:size(I,1)+1:size(I,1)*size(I,1))= normalisation;

imdl.inv_solve_abs_core.calc_solution_error = 0;

imdl.inv_solve_abs_core.verbose = 10;
%%imdl.inv_solve_abs_core.normalisation= I;
%%imdl.inv_solve_abs_core.homogeneization= 1;
% imdl.inv_solve_abs_core.fixed_background= 1; % the default now in _core
imdl.inv_solve_abs_core.line_search_args.perturb= [0 5*logspace(-7,-4,5)];
%imdl.inv_solve_abs_core.max_iterations= 10;
%imdl.inv_solve_abs_core.plot_line_optimize = 1;

imdl.inv_solve_abs_core.elem_output = 'log10_resistivity';
imgr= inv_solve_abs_core(imdl, dd);

% save the result so we don't have to wait forever if we want to look at the result later
%save('inv_solve_abs_core_rec2.mat', 'imgr');

imgGNd= imgr;
%imgGNd.fwd_model.coarse2fine= cmdl.coarse2fine;
% removal of the background elem_data is now handled in the solver
% conversion to output elem_data is now handled in the solver
%imgGNd.elem_data= log10(imgGNd.res_data(1:end-1));
%imgGNd.calc_colours.clim= 1.5;
%imgGNd.calc_colours.ref_level= 1.5;

elec_posn= zeros(length(fmdl.electrode),3);
for i=1:length(fmdl.electrode)
    elec_posn(i,:)= mean(fmdl.nodes(fmdl.electrode(1,i).nodes,:),1);
end

clf; show_fem(imgGNd,1);
hold on; plot(elec_posn(:,1),elec_posn(:,3),'k*');
axis tight; ylim([-100 0.5])
xlabel('X (m)','fontsize',20,'fontname','Times')
ylabel('Z (m)','fontsize',20,'fontname','Times')
set(gca,'fontsize',20,'fontname','Times');

img = mk_image( imdl );
img.elem_data= 1./(10.^imgr.elem_data);
vCG= fwd_solve(img); vCG = vCG.meas;

I = 1; % TODO FIXME -> I is diag(1./vh) the conversion to apparent resistivity
% TODO these plots are useful, get them built into the solver!
clf; plot(I*(dd.meas-vCG)); title('data misfit');
clf; hist(abs(I*(dd.meas-vCG)),50); title('|data misfit|, histogram'); xlabel('|misfit|'); ylabel('count');
clf; show_pseudosection( fmdl, I*dd.meas); title('measurement data');
clf; show_pseudosection( fmdl, I*vCG); title('reconstruction data');
clf; show_pseudosection( fmdl, (vCG-dd.meas)./dd.meas*100); title('data misfit');
