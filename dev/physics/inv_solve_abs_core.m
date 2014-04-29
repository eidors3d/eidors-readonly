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
% Parameters (inv_model.parameters.*):
%   verbose (show progress)                (default 4)
%      0: quiet
%    >=1: print iteration count
%    >=2: print details as the algorithm progresses
%    >=3: plot residuals versus iteration count
%    >=4: plot line search per iteration
%   plot_residuals                         (default 0)
%    plot residuals without verbose output
%   fwd_solutions                          (default 0)
%    0: ignore
%    1: count fwd_solve(), generally the most
%       computationally expensive component of
%       the iterations
%   residual_func =             (default @GN_residual)
%    NOTE: @meas_residual exists to maintain
%    compatibility with some older code
%   beta_func =                        (default: none)
%   max_iterations                        (default 10)
%   ntol (estimate of machine precision) (default eps)
%   tol (stop iter if r_k < tol)           (default 0)
%   dtol                              (default -0.01%)
%    stop iter if (r_k - r_{k-1}) < dtol AND
%                 k >= dtol_iter
%   dtol_iter                              (default 0)
%    apply dtol stopping criteria if k >= dtol_iter
%   min_value                           (default -inf)
%   max_value                           (default +inf)
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
%    background value is required
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
%   elem_working              (default 'conductivity')
%   elem_output               (default 'conductivity')
%    The working and output units for 'elem_data'.
%    Valid types are 'conductivity' and 'resistivity'
%    as plain units or with the prefix 'log_' or
%    'log10_'. Conversions are handled internally.
%    Scaling factors are applied to the Jacobian
%    (calculated in units of 'conductivity') as
%    appropriate, see calc_jacobian_scaling_func.
%    If elem_working == elem_output, then no
%    conversions take place.
%   meas_input                     (default 'voltage')
%   meas_working                   (default 'voltage')
%    Similarly to elem_working/output, conversion
%    between 'voltage' and 'apparent_resistivity' and
%    their log/log10 varients are handled internally.
%    If meas_input == meas_working no conversions take
%    place. The normalization factor 'N' is calculated
%    if 'apparent_resistivity' is used.
%   calc_jacobian_scaling_func           (default ...)
%    See elem_working for defaults. This is used to
%    apply the chain rule to a Jacobian calculated for
%    conductivity. If the Jacobian in use is not the
%    default, a scaling function is required.
%   normalize_data_func
%    (default  @apparent_resistivity_factor)
%    dv = N*(data-data0) where N is provided by this
%    function
%
%   Signature for residual_func
%    r = f(dv, de, W, hp2, RtR)
%   where
%    r   - the residual
%    dv  - change in voltage
%    de  - change in image elements
%    W   - measurement inverse covarience matrix
%    hp2 - hyperparameter squared, see CALC_HYPERPARAMETER
%    RtR - regularization matrix squared
%
%   Signature for line_optimize_func
%    [alpha, img, dv, opt] = f(img, sx, data0, img0, N, W, hp2, RtR, dv, opt)
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
%    RtR   - regularization matrix squared
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
%   Signature for  normalize_data_func
%    N = f(fmdl)
%   where
%    N    - the normalization factor
%    fmdl - the forward model, used to estimate the N
%    when voltage to apparent_resistivity is required
%    This function can be repurposed for scalingi
%    measurement data in other ways.
%
%   Signature for beta_func
%   (zero if not supplied, for use in Conjugate Gradient)
%    beta = f(dx_k, dx_{k-1}, sx_{k-1})
%   where
%    dx_k     - the current descent direction
%    dx_{k-1} - the previous descent direction
%    sx_{k-1} - the previous search direction (after applying beta)
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
% (C) 2010-2014 Andy Adler & Bartłomiej Grychtol, Nolwenn Lespare, Alistair Boyle.
% License: GPL version 2 or version 3

% $Id$

%AB->BG: The key item to discuss is my attempt to attack the
% parameterization issue. This is by no means a generic
% solution but it does present a very clean interface to the
% user. I recognize this implementation does not quite fit
% with your view of how parameterization should exist in the
% rest of EIDORS. What are your thougths? I am open to any
% solution/refactoring that would maintain the same
% functionality/interface here or provide a clean
% equivalent.

%AB->BG: It is my intent to not expose this function
% directly but to have a inv_solve_abs_GN() wrapper.
% I have constructed this with the intent that it is trivial
% to implement a CG solver using this framework. Hopefully,
% I'll get to that shortly. This would result in a pair of
% inv_solve_abs_CG and inv_solve_abs_CG_log wrapper function.
% This file is a bit of a monster but it is removing a bit
% of cut and paste that was expanding rapidly.
% Do you have thoughts about this approach?

%AB->BG: The UNIT_TEST is working. It does not seem to be
% nearly challenging enough to show much of the odd behaviour
% I spent two months debugging in the fall. I am not sure
% what to do about this. It seems that building a truely
% robust framework here requires some tough UNIT_TESTs so
% that the infrastructure will not be broken as it gets
% tweaked. Do you have thoughts or ideas about how to go
% about this or how far you think UNIT_TESTs should go?
% As an example, I believe they need to be self-checking to
% be truly useful in the long run but that is pretty tough
% to construct in a bullet proof manner.

%AB->BG: There are still TODO comments in the code and
% further clean up that could/should be done. I am holding
% off on this until we can discuss and/or resolve some of the
% outstanding external issues, as above and also some of the
% TODO comments bring issues to mind that we have discussed
% before such as the mk_image()'s handling of c2f.

%AB->BG: elem_movement_init (not documented above) was my
% only concession to electrode movement: we need an electrode
% movement 'background' constructor. For the rest of the
% code, the movement is handled as 'just another odd
% parameter that doesn't map to conductivity'. A clean way
% to handle this as part of the 'parameterization' concept
% is important.


%--------------------------
% UNIT_TEST?
if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST') && (nargin == 1); img = do_unit_test; return; end
if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST') && (nargin == 2); img = do_unit_test(data0); return; end

%--------------------------
opt = parse_options(inv_model);
inv_model = stupid_model_translation(opt, inv_model); % TODO rm -- transitional function to DIE
%if opt.do_starting_estimate
%    img = initial_estimate( inv_model, data0 ); % TODO
%%%    AB->NL this is Nolwenn's homogeneous estimate...
%%%    calc_background_resistivity is my version of this code
%%%    that is working for my data set
%else
[inv_model, opt] = append_c2f_background(inv_model, opt);
img = mk_image( inv_model );
%img = calc_jacobian_bkgnd( inv_model );
% TODO does calc_jacobian_bkgnd ignore 'physics' right now.. that might screw things up pretty good!
% img = physics_param_mapper(img); % copy data from whatever 'physics' to img.params
% mk_image doesn't handle the c2f
% TODO move this from here to mk_image so we get an img with the correct number of elem_data
if isfield(inv_model, 'fwd_model') && ...
   isfield(inv_model.fwd_model, 'coarse2fine') && ...
   length(img.(opt.elem_working).elem_data) ~= size(inv_model.fwd_model.coarse2fine,2)

   % TODO replace with fix_c2f_calc_jacobian_backgnd
   c2f = inv_model.fwd_model.coarse2fine;
   %img.elem_data = c2f \ img.elem_data; % --> rank deficient???
   bg = mean(img.(opt.elem_working).elem_data);
   img.(opt.elem_working).elem_data = ones(size(c2f,2),1)*bg;
   if opt.verbose > 1
      bg_tmp = map_img(bg, opt.elem_working, 'resistivity');
      fprintf('  c2f: correcting mk_image elem_data size %d -> %d (av %0.1f Ohm.m)\n', size(c2f), bg_tmp);
      disp('    TODO this fix should be moved to mk_image()');
   end
end
% transfer parameter into img
if isfield(opt, 'elem_movement_init')
   img.elem_movement_init = opt.elem_movement_init;
end

% precalculate some of our matrices if required
hp2 = init_hp(inv_model, opt);
W  = init_meas_icov(inv_model, opt);
N = init_normalization(inv_model.fwd_model, opt);

% map data and measurements to working types
%  convert elem_data
% img = map_img(img, opt.elem_working); %DIE

data0 = stupid_data_translation(opt, data0, N); % TODO rm -- transitional function to DIE

% now get on with
img0 = physics_param_mapper(img);
RtR = 0; k = 0; dv = []; de = []; sx = 0; r = 0; stop = 0; % general init
residuals = zeros(opt.max_iterations,3); fig_r = []; % for residuals plots
if opt.verbose > 1
   fprintf('  iteration start up\n')
end
dxp = 0;
while 1
  % update RtR, if required (depends on prior)
%   img = physics_param_mapper(img);
  RtR = update_RtR(RtR, inv_model, k, img, opt);

  % update change in element data from the prior de and
  % the measurement error dv
  [dv, opt] = update_dv(dv, img, data0, N, opt);
  de = update_de(de, img, img0, opt);

  % now find the residual, quit if we're done
  [stop, k, r, fig_r] = update_residual(dv, de, W, hp2, RtR, k, r, fig_r, opt);
  if stop
     break;
  end
  if opt.verbose > 1
     fprintf('  iteration %d\n', k)
  end

  % calculate the Jacobian
  J = update_jacobian(img, N, opt);

  % determine the next search direction sx
  %  dx is specific to the algorithm, generally "downhill"
  dx = update_dx(J, W, hp2, RtR, dv, de, opt);
  % choose beta, beta=0 unless doing Conjugate Gradient
  beta = update_beta(dx, dxp, sx, opt);
  % sx_k = dx_k + beta * sx_{k-1}
  sx = update_sx(dx, beta, sx, opt);
  dxp = dx; % saved for next iteration if using beta

  % line search for alpha, leaving the final selection as img
  [alpha, img, dv, opt] = update_alpha(img, sx, data0, img0, N, W, hp2, RtR, dv, opt);
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
% convert data for output
img = physics_data_mapper(img);
img = rmfield(img, opt.elem_working);
img = map_img(img, opt.elem_output);
img.(opt.elem_output).elem_data = [];
img = physics_data_mapper(img,1); % move data from img.elem_data to whatever 'physics'

function W = init_meas_icov(inv_model, opt)
   W = 1;
   if opt.calc_meas_icov
      if opt.verbose > 1
         disp('  calc measurement inverse covariance W');
      end
      W   = calc_meas_icov( inv_model );
   end

function hp2 = init_hp(inv_model, opt)
   hp2 = 0;
   if opt.calc_hyperparameter
      if opt.verbose > 1
         disp('  calc regularization hyperparameter(s)');
      end
      hp2  = calc_hyperparameter( inv_model ).^2;
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
function [stop, k, r, fig_r] = update_residual(dv, de, W, hp2, RtR, k, r, fig_r, opt)
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
  r_k = feval(opt.residual_func, dv, de, W, hp2, RtR);
  % save residual for next iteration
  r(k,1) = r_k;

  % now do something with that information
  if opt.verbose > 1
     if k == 1
        fprintf('    calc residual, r=%0.3g\n', r_k);
     else
        fprintf('    calc residual\n');
        fprintf('      r =%0.3g\n', r_k);
        dr = (r_k - r_km1);
        fprintf('      dr=%0.3g (%0.3g%%)\n', dr, dr/r_km1*100);
     end
  end
  if opt.plot_residuals
     %         optimization_criteria, data misfit, roughness
     r(k,2:3) = [(dv'*dv)/2 (de'*de)/2];
     if k > 1
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
  end

  % evaluate stopping criteria
  if k > opt.max_iterations
     stop = 1;
  end
  % TODO return 'measurement residual' & 'roughness' for progress plot, as well
  if r_k < opt.tol + opt.ntol
     if opt.verbose > 1
        fprintf('  terminated at iteration %d\n',k);
        fprintf('    residual tolerance (%0.3g) achieved\n', opt.tol + opt.ntol);
     end
     stop = 1;
  end
  if (k > opt.dtol_iter) && ((r_k - r_km1)/r_km1 > opt.dtol + 2*opt.ntol)
     if opt.verbose > 1
        fprintf('  terminated at iteration %d (iterations not improving)\n', k);
        fprintf('    residual slope tolerance (%0.3g) exceeded\n', opt.dtol + 2*opt.ntol);
     end
     stop = 1;
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
      beta = feval(opt.beta_func, dx_k, dx_km1, sx_km1);
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
         fprintf( '      direction change ||ddx||=%0.3g\n', norm(sx/nsx-sx_km1/nsxk));
      end
   end

function RtR = update_RtR(RtR, inv_model, k, img, opt)
   % TODO sometimes (with Noser?) this requires the Jacobian, could this be done more efficiently?
   % add a test function to determine if img.elem_data affects RtR, skip this if independant
   % TODO we could detect in the opt_parsing whether the calc_RtR_prior depends on 'x' and skip this if no
   if opt.calc_RtR_prior
      if opt.verbose > 1
         try RtR_str = func2str(inv_model.RtR_prior);
         catch
            try RtR_str = inv_model.RtR_prior;
            catch RtR_str = 'unknown';
            end
         end
         fprintf('    calc regularization RtR (%s)\n', RtR_str);
      end
      RtR = calc_RtR_prior_wrapper(inv_model, img, opt);
   else
      RtR = 0;
   end

function RtR = calc_RtR_prior_wrapper(inv_model, img, opt)
%    img = map_img(img, opt.elem_working);
   img = physics_param_mapper(img);
   inv_model.jacobian_backgnd = img;
   RtR = calc_RtR_prior( inv_model );
   if size(RtR,1) < length(img.params)
     ne = length(img.params) - size(RtR,1);
     % we are correcting for the added background element
     % if there is movement, don't add it there.
     if isfield(img, 'current_physics') && ...
        any(strcmp(img.current_physics, 'movement'))
        % grab the first set of non-movement elem_data
        i = find(~strcmp(img.current_physics, 'movement'));
        ps = img.physics_sel{i(1)};
        % insert an extra element at the end of the non-movement
        % regularization so the length is right, we assume Tikhonov
        RtR(ps(end)+1+ne:end+ne, ps(end)+1+ne:end+ne) = RtR(ps(end)+1:end, ps(end)+1:end);
        RtR(ps(end)+1:ps(end)+ne, ps(end)+1:ps(end)+ne) = RtR(ps(1),ps(1));
        if opt.verbose > 1
           fprintf('    c2f: adjusting RtR by appending %d rows/cols to non-movement\n', ne);
        end
     else
        RtR(end+1:end+ne, end+1:end+ne) = RtR(1,1);
        if opt.verbose > 1
           fprintf('    c2f: adjusting RtR by appending %d rows/cols\n', ne);
        end
     end
     disp('      TODO move this fix, or something like it to calc_RtR_prior -- this fix is a quick HACK to get things to run...');
   end

function J = update_jacobian(img, N, opt)   
   if(opt.verbose > 1)
      try J_str = func2str(img.fwd_model.jacobian);
      catch J_str = img.fwd_model.jacobian;
      end
      fprintf('    calc Jacobian J(x) (%s,', J_str);
   end
%DIE    if donew
%DIE      img.fwd_model.measured_quantity = 'apparent_resistivity';
     J = calc_jacobian(img);
%DIE    else
%DIE      img = map_img(img, 'conductivity');   
%DIE      % scaling if we are working in something other than direct conductivity
%DIE      S = feval(opt.calc_jacobian_scaling_func, img.elem_data); % chain rule
%DIE      % finalize the jacobian
%DIE      % Note that if a normalization (i.e. apparent_resistivity) has been applied
%DIE      % to the measurements, it needs to be applied to the Jacobian as well!
%DIE      J = N * calc_jacobian( img ) * S;
%DIE    end
   if opt.verbose > 1
      fprintf(' %d DoF, %d meas, %s)\n', size(J,2)-length(opt.elem_fixed), size(J,1), func2str(opt.calc_jacobian_scaling_func));
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


function [alpha, img, dv, opt] = update_alpha(img, sx, data0, img0, N, W, hp2, RtR, dv, opt)
  if(opt.verbose > 1)
     disp('    line search');
  end
  tmp = physics_param_mapper(img);
  % some sanity checks before we feed this information to the line search
  err_if_inf_or_nan(sx, 'sx (pre-line search)');
  err_if_inf_or_nan(tmp.params, 'img.params (pre-line search)');

  if any(size(tmp.params) ~= size(sx))
     error(sprintf('mismatch on params[%d,%d] vs. sx[%d,%d] vector sizes, check c2f_background_fixed',size(img.elem_data), size(sx)));
  end
  [alpha, img, dv, opt] = feval(opt.line_search_func, img, sx, data0, img0, N, W, hp2, RtR, dv, opt);
  if(opt.verbose > 1)
     fprintf('      selected alpha=%0.3g\n', alpha);
  end

function err_if_inf_or_nan(x, str);
  if any(isnan(x) | isinf(x))
      error(sprintf('bad %s (%d NaN, %d Inf of %d)', ...
                    str, ...
                    length(find(isnan(x))), ...
                    length(find(isinf(x))), ...
                    length(x)));
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
%    img0 = map_img(img0, opt.elem_working);
%    img  = map_img(img,  opt.elem_working);
   img = physics_param_mapper(img);
   err_if_inf_or_nan(img0.params, 'de img0');
   err_if_inf_or_nan(img.params,  'de img');
   % probably not the most robust check for whether this is the first update
   % but this ensures that we get exactly zero for the first iteration and not
   % a set of values that has numeric floating point errors that are nearly zero
   if isempty(de) % first iteration
      % data hasn't changed yet!
      de = zeros(size(img0.params));
   else
      de = img0.params - img.params;
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
   [dv, opt] = update_dv_core(img, data0, N, opt);

% also used by the line search as opt.line_search_dv_func
function [dv, opt] = update_dv_core(img, data0, N, opt)
%    img = map_img(img, 'conductivity'); %DIE Why can't I move this?
   img.fwd_model.measured_quantity = opt.meas_working;
   data = fwd_solve(img);
   dv = calc_difference_data(data, data0, img.fwd_model);
   err_if_inf_or_nan(dv, 'dv out');

function do = donew; do =1;

function show_fem_iter(k, img, inv_model, opt)
  if opt.verbose > 1
     disp('    show_fem()');
  end
  out = opt.elem_output;
  if iscell(out)
     out = out{1};
  end
  img = map_img(img, out);
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
  figure; show_fem(img, 1);
  str = strrep(out, '_', ' ');
  title(sprintf('iter=%d, %s',k, str));

% TODO confirm that GN line_search_onm2 is using this residual calculation (preferably, directly)
function residual = GN_residual(dv, de, W, hp2, RtR)
%   [size(dv); size(W); size(de); size(hp2RtR)]
   % we operate on whatever the iterations operate on (log data, resistance, etc) + perturb(i)*dx
   hp2RtR = hp2*RtR;
   residual = 0.5*( dv'*W*dv + de'*hp2RtR*de);

function residual = meas_residual(dv, de, W, hp2, RtR)
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

function imdl = stupid_model_translation(opt, imdl)
  imdl.fwd_model.measured_quantity = opt.meas_working;
  if donew
    if isfield(imdl.jacobian_bkgnd, 'value') % TODO: do this smarter
     bkg = mk_image(imdl.fwd_model,imdl.jacobian_bkgnd.value);
     bkg = map_img(bkg,opt.elem_working);
     bkg.(opt.elem_working).elem_data = bkg.elem_data;
     bkg.current_physics = [];
     bkg = rmfield(bkg,{'fwd_model','elem_data'});
     imdl.jacobian_bkgnd = bkg;
    end
  else
    % we're okay thanks...
  end

function data0 = stupid_data_translation(opt, data0, N)
  data0 = map_meas(data0, opt.meas_working, N);

  if donew
  else
    %  nothing
  end

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
      opt.residual_func = @GN_residual; % r = f(dv, de, W, hp2, RtR)
      % NOTE: the meas_residual function exists to maintain
      % compatibility with Nolwenn's code, the GN_residual
      % is a better choice
      %opt.residual_func = @meas_residual; % r = f(dv, de, W, hp2, RtR)
   end

   % calculation of update components
   if ~isfield(opt, 'update_func')
      opt.update_func = @GN_update; % dx = f(J, W, hp2, RtR, dv, de)
   end
   % figure out if things need to be calculated
   if ~isfield(opt, 'calc_meas_icov') % derivative of the objective function
      opt.calc_meas_icov = 0; % W
   end
   if ~isfield(opt, 'calc_RtR_prior') % derivative of the objective function
      opt.calc_RtR_prior = 0; % RtR
   end
   if ~isfield(opt, 'calc_hyperparameter') % derivative of the objective function
      opt.calc_hyperparameter = 0; % hp2
   end
%   try
      if opt.verbose > 1
         fprintf('    examining function %s(...) for required arguments\n', func2str(opt.update_func));
      end
      % ensure that necessary components are calculated
      % opt.update_func: dx = f(J, W, hp2, RtR, dv, de)
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
      disp('    inv_model.parameters.plot_residuals=0');
   end

   % line search
   if ~isfield(opt,'line_search_func')
      % [alpha, img, dv, opt] = f(img, sx, data0, img0, N, W, hp2, RtR, dv, opt);
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

   % meas_select already handles selecting from the valid measurements
   % we want the same for the elem_data, so we only work on modifying the legal values
   % Note that c2f_background's elements are added to this list if opt.c2f_background_fixed == 1
   if ~isfield(opt, 'elem_fixed')
      opt.elem_fixed = [];
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
      % TODO error out if elem_working is not onductivity and the jacobian function is not the default... we can't guess correctly then and we'll get funky/had-to-debug behaviour
      if ~strcmp(opt.elem_working, 'conductivity') && ...
         (~isfield(imdl, 'fwd_model') || ...
          ~isfield(imdl.fwd_model, 'jacobian') || ...
          ~strcmp(imdl.fwd_model.jacobian, 'eidors_default'))
         error('can not guess at inv_model.parameters.calc_jacobian_scaling_func, one must be provided');
      end
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
      if any(strcmp({opt.meas_input, opt.meas_working}, 'apparent_resistivity'))
         opt.normalize_data_func = @apparent_resistivity_factor; % N = f(fmdl)
      end
   end
   if isfield(opt, 'normalize_data_func')
      opt.normalize_data = 1;
   else
      opt.normalize_data = 0;
   end

function check_matrix_sizes(J, W, hp2, RtR, dv, de, opt)
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
      any(size(RtR) ~= [ne ne])
      error('RtR size (%d rows, %d cols) is incorrect (%d rows, %d cols)', size(RtR), ne, ne);
   end

function dx = update_dx(J, W, hp2, RtR, dv, de, opt)
   if(opt.verbose > 1)
      fprintf( '    calc step size dx');
   end

   % don't penalize for fixed elements
   de(opt.elem_fixed) = 0;

   % TODO move this outside the inner loop of the iterations, it only needs to be done once
   check_matrix_sizes(J, W, hp2, RtR, dv, de, opt)

   % do the update step direction calculation
   dx = feval(opt.update_func, J, W, hp2, RtR, dv, de);
   % ignore any fixed value elements
   dx(opt.elem_fixed) = 0;

   if(opt.verbose > 1)
      fprintf(', ||dx||=%0.3g\n', norm(dx));
   end

function dx = GN_update(J, W, hp2, RtR, dv, de)
   hp2RtR = hp2*RtR;
   % the actual update
   dx = (J'*W*J + hp2RtR)\(J'*dv + hp2RtR*de);

% for each argument, returns 1 if the function depends on it, 0 otherwise
% 'zero' arguments do not need to be calculated since they don't get used
function args = function_depends_upon(func, argn)
   % build function call
   str = 'feval(func,';
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

% this function always returns zero
function out = ret0_func(arg1, arg2);
   out = 0;

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
         opt.elem_fixed(end+1) = nc+1;
      end
    end

function [img, opt] = strip_c2f_background(img, opt, indent)
    if nargin < 3
       indent = '';
    end
    % nothing to do?
    if opt.c2f_background <= 0
      return;
    end

    % if there are multiple 'physics', we assume its the first
    % TODO -- this isn't a great assumption but it'll work for now,
    %         we should add a better (more general) mechanism
    in = opt.elem_working;
    out = opt.elem_output;
    if iscell(in)
       in = in{1};
    end
    if iscell(out)
       out = out{1};
    end

    % go about cleaning up the background
    e = opt.c2f_background;
    % take backgtround elements and convert to output 'physics' (resistivity, etc)
    img = physics_param_mapper(img);
    bg = map_img(img.params(e), in, out);

    img.elem_data_background = bg;
    % remove elements from elem_data & c2f
    img.params(e) = [];
    img.fwd_model.coarse2fine(:,e) = [];
    % remove our element from the lists
    opt.c2f_background = 0;
    ri = find(opt.elem_fixed == e);
    opt.elem_fixed(ri) = [];
    if isfield(img, 'physics_sel')
       for i = 1:length(img.physics_sel)
          t = img.physics_sel{i};
          ti = find(t == e);
          t(ti) = []; % rm 'e' from the list of physics_sel
          ti = find(t > e);
          t(ti) = t(ti)-1; % down-count element indices greater than our deleted one
          img.physics_sel{i} = t;
       end
    end

    % show what we got for a background value
    if(opt.verbose > 1)
       bg = img.elem_data_background;
       bg = map_img(bg, in, 'resistivity');
       fprintf('%s  background conductivity: %0.1f Ohm.m\n', indent, bg);
    end

    img = physics_param_mapper(img,1);
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

function pass = do_unit_test(solver)
if nargin == 0
  solver = 'inv_solve_abs_core';
end
pass = 1;
pass = pass & do_unit_test_sub;
pass = pass & do_unit_test_rec1(solver);
%pass = pass & do_unit_test_rec2; % TODO this unit test is very, very slow... what can we do to speed it up... looks like the perturbations get kinda borked when using the line_search_onm2
if pass
   disp('TEST: overall PASS');
else
   disp('TEST: overall FAIL');
end

% test sub-functions
% jacobian scalings
function pass = do_unit_test_sub
pass = 1;
d = 1;
while d ~= 1 & d ~= 0
  d = rand(1);
end

disp('TEST: Jacobian scaling');
d = [d d]';
if any(ret1_func(d) ~= 1)
   fprintf('TEST: FAIL for Jacobian scaling (%s)\n', 'conductivity');
   pass = 0;
end
if any(dx_dlogx(d) ~= diag(d))
   fprintf('TEST: FAIL for Jacobian scaling (%s)\n', 'log_conductivity');
   pass = 0;
end
if any(dx_dlog10x(d) ~= diag(d)*log(10))
   fprintf('TEST: FAIL for Jacobian scaling (%s)\n', 'log10_conductivity');
   pass = 0;
end
if any(dx_dy(d) ~= diag(-d.^2))
   fprintf('TEST: FAIL for Jacobian scaling (%s)\n', 'resistivity');
   pass = 0;
end
if any(dx_dlogy(d) ~= diag(-d))
   fprintf('TEST: FAIL for Jacobian scaling (%s)\n', 'log_resistivity');
   pass = 0;
end
if any(dx_dlog10y(d) ~= diag(-d)/log(10))
   fprintf('TEST: FAIL for Jacobian scaling (%s)\n', 'log10_resistivity');
   pass = 0;
end

% a couple easy reconstructions
% check c2f, apparent_resistivity, log_conductivity, verbosity don't error out
function pass = do_unit_test_rec1(solver)
pass = 1;
% -------------
% ADAPTED FROM
% Create simulation data $Id: basic_iterative01.m 3829 2013-04-13 14:21:30Z bgrychtol $
%  http://eidors3d.sourceforge.net/tutorial/adv_image_reconst/basic_iterative.shtml
% 3D Model
imdl= mk_common_model('c2t4',16); % 576 elements
imdl.solve = solver;
imdl.reconst_type = 'absolute';
imdl.parameters.elem_working = 'log_conductivity';
imdl.parameters.meas_working = 'apparent_resistivity';
imdl.inv_solve.calc_solution_error = 0;
imdl.parameters.verbose = 0;
%show_fem(imdl.fwd_model);
imgsrc= mk_image( imdl.fwd_model,1, 'conductivity');
% set homogeneous conductivity and simulate
vh=fwd_solve(imgsrc);
imgsrc = physics_data_mapper(imgsrc); %DIE: so we can work on elem_data
% set inhomogeneous conductivity and simulate
ctrs= interp_mesh(imdl.fwd_model);
x= ctrs(:,1); y= ctrs(:,2);
r1=sqrt((x+5).^2 + (y+5).^2); r2 = sqrt((x-85).^2 + (y-65).^2);
imgsrc.elem_data(r1<50)= 0.05;
imgsrc.elem_data(r2<30)= 100;
imgp = map_img(imgsrc, 'log10_conductivity');
imgsrc = physics_data_mapper(imgsrc,1); % DIE: revert
imgp = rmfield(imgp,'conductivity');
imgp.log10_conductivity.elem_data = []; %DIE: physics_data_mapper needs it
imgp = physics_data_mapper(imgp,1); % revert

hh=gcf; subplot(221); show_fem(imgp,1); axis tight; title('synthetic data, logC');
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
%imdl.parameters.meas_working = 'apparent_resistivity';
img2= inv_solve(imdl, vi);
figure(hh); subplot(223); show_fem(img2,1); axis tight; title('#2 verbosity=0');
img1 = physics_data_mapper(img1);
img2 = physics_data_mapper(img2);
max_err = max(abs((img1.elem_data - img2.elem_data)./(img1.elem_data)));
if max_err > 0.05
  fprintf('TEST:  img1 != img2 --> FAIL %g %%\n', max_err);
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
e1 = c2f \ img1.elem_data; % noisy and unstable... but its a crude check
img3 = physics_data_mapper(img3);
e3 = img3.elem_data;
err = abs((e1 - e3) ./ e1);
err(abs(e1) < 20) = 0;
err_thres = 0.40;
if any(err > err_thres) % maximum 15% error
  ni = find(err > err_thres);
  fprintf('TEST:  img1 != img3 --> FAIL max(err) = %0.2e on %d elements (thres=%0.2e)\n', ...
          max(err(ni)), length(ni), err_thres);
  pass = 0;
else
  disp('TEST:  img1 == img3 --> PASS');
end

%imdl.parameters.verbose = 1000;
imdl.parameters.elem_output = 'log10_resistivity'; % resistivity output works
img4= inv_solve(imdl, vi);
figure(hh); subplot(224); show_fem(img4,1); axis tight; title('#4 c2f + log10 resistivity out');
% check
img4 = physics_data_mapper(img4);
e4 = 1./(10.^img4.elem_data);
err = abs((e1 - e4) ./ e1);
err(abs(e1) < 20) = 0;
err_thres = 0.40;
if any(err > err_thres) % maximum 15% error
  ni = find(err > err_thres);
  fprintf('TEST:  img1 != img4 --> FAIL max(err) = %0.2e on %d elements (thres=%0.2e)\n', ...
          max(err(ni)), length(ni), err_thres);
  pass = 0;
else
  disp('TEST:  img1 == img4 --> PASS');
end


function pass = do_unit_test_rec2
disp('TEST: reconstruct a discontinuity');
pass = 1; % TODO fail criteria?
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
figure; show_fem(img); title('model');

% img2= mk_image(fmdl,img.elem_data);
% figure; show_fem(img2);

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
imdl.solve = @inv_solve_abs_core;
imdl.reconst_type = 'absolute';
imdl.hyperparameter.value = 1e2; % was 0.1
imdl.jacobian_bkgnd.value = 1;

imdl.parameters.elem_working = 'log_conductivity';
imdl.parameters.meas_working = 'apparent_resistivity';
imdl.parameters.dtol_iter = 4; % default 1 -> start checking on the first iter
imdl.parameters.max_iterations = 20; % default 10

% the conversion to apparaent resistivity is now handled inside the solver
%%img1= mk_image(fmdl,1);
%%vh1= fwd_solve(img1);
%%normalisation= 1./vh1.meas;
%%I= speye(length(normalisation));
%%I(1:size(I,1)+1:size(I,1)*size(I,1))= normalisation;

imdl.inv_solve.calc_solution_error = 0;

imdl.parameters.verbose = 10;
%%imdl.parameters.normalisation= I;
%%imdl.parameters.homogeneization= 1;
% imdl.parameters.fixed_background= 1; % the default now in _core
imdl.parameters.line_search_args.perturb= [0 5*logspace(-7,-4,5)];
%imdl.parameters.max_iterations= 10;
%imdl.parameters.plot_line_optimize = 1;

imdl.parameters.elem_output = 'log10_resistivity';
imgr= inv_solve(imdl, dd);

% save the result so we don't have to wait forever if we want to look at the result later
save('inv_solve_abs_core_rec2.mat', 'imgr');

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

figure; show_fem(imgGNd,1);
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
figure; plot(I*(dd.meas-vCG)); title('data misfit');
figure; hist(abs(I*(dd.meas-vCG)),50); title('|data misfit|, histogram'); xlabel('|misfit|'); ylabel('count');

figure; show_pseudosection( fmdl, I*dd.meas); title('measurement data');
figure; show_pseudosection( fmdl, I*vCG); title('reconstruction data');
figure; show_pseudosection( fmdl, (vCG-dd.meas)./dd.meas*100); title('data misfit');
