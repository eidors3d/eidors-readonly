function  [alpha, img, dv, opt] = line_search_onm2(imgk, dx, data1, img1, N, W, hp2RtR, dv0, opt, retry, pf_max)
% function  [alpha, img, dv, opt] = line_search_onm2(imgk, dx, data1, img1, N, W, hp2RtR, dv0, opt)
% line search function with a fitted polynomial of O(n-2) where n is the number of perturbations
% (C) 2013 Alistair Boyle
% License: GPL version 2 or version 3

if nargin < 11
  retry = 0;
end
perturb= sort(opt.line_search_args.perturb);
if nargin < 12
  pf_max = length(perturb)-2;
end

if opt.verbose > 1
   fprintf('     ');
end
% fwd_solve is the most expensive part generally, count how many we do
if ~isfield(opt, 'fwd_solutions')
   opt.fwd_solutions = 0;
end
x = imgk.elem_data;

if(perturb(1) ~= 0)
  error('line_search_o2() expects first perturbation (inv_model.parameters.line_search.perturb) to be alpha=0');
end
% Compute the forward model for each perturbation step
img = imgk;
% mlist is our search result for each alpha value, perturb(i)
%  -- NaN: initiailized but not calculated
%  -- -Inf: should not occur we have code that converts calculated NaNs and -Inf to +Inf
%  -- +Inf: calculated value was bad, ignore it
mlist= ones(size(perturb))*NaN; % init
idx_list = [ 1 length(perturb):-1:2 ];
for i = idx_list;
    if opt.verbose > 1
       fprintf(' [%d]=%0.3g', i, perturb(i));
    end
    if (i == 1) && (~isempty(dv0))
      % don't bother simulating when alpha=0 (we already have the measurements)
      dv = dv0; % @ alpha=0 from the previous line search iteration
    else
      % fwd_solve and then calculate measurement error (dv)
      img.elem_data = x + perturb(i)*dx;
      [dv, opt] = feval(opt.line_search_dv_func, img, data1, N, opt);
    end
    % build de, the change in elem_data from the initial guess
    de = feval(opt.line_search_de_func, img, img1, opt);
    % we only calculate a new residual if the input data is "sane"
    if any(isnan(dv) | isinf(dv))
       warning(sprintf('%d of %d elements in dv are NaN or Inf', ...
                       length(dv), ...
                       length(find(isnan(dv) | isinf(dv)))));
       mlist(i) = +Inf;
    elseif any(isnan(de) | isinf(de))
       warning(sprintf('%d of %d elements in de are NaN or Inf', ...
                       length(de), ...
                       length(find(isnan(de) | isinf(de)))));
       mlist(i) = +Inf;
    else
       % calculate the residual
       mlist(i) = feval(opt.residual_func, dv, de, W, hp2RtR);
       if any(isnan(mlist(i)) | isinf(mlist(i)))
          mlist(i) = +Inf; % NaN or Inf are converted to Inf, since we use NaN to indicate initialized but not calculated
       end
    end
end
if opt.verbose > 1
   fprintf('\n');
   fprintf('      fitting data\n      ');
   fprintf('  %0.3g',mlist);
   fprintf('\n');
end
% drop bad values
if isinf(mlist) % NaN's from any calculation were converted to Inf's earlier
   warning('encoutered NaN or +-Inf residuals, something has gone wrong in the line search, converting to large numbers and carrying on');
end

% TODO looks like this was quiting when there are still good choices left
%if all((mlist/mlist(1)-1) < 1e-4) % < 0.01% change
%   % TODO maybe we need to search *larger* perturbations here... for now we just short circuit the repeated retries at the end, when we are not improving
%   if opt.verbose > 1
%      fprintf('      stopping line search: no further improvements observed\n');
%   end
%   img = imgk;
%   alpha = 0;
%   dv = dv0;
%   return;
%end

% For our poly fit, we drop all NaN and Inf values
goodi = find((~isnan(mlist)) & (~isinf(mlist)));
alpha=perturb(end);
meas_err = +Inf; % make sure we grab the min(mlist) if we're not doing a polyfit
if length(goodi) > 2
  % Select the best fitting step, we scale and
  p_rng = range(perturb(goodi)); % p_min = 0
  pfx = log10(perturb(goodi)/p_rng);
  pfx(1) = -100; % log10(0) = -Inf --> -1e100 so that it's finite
  pf= polyfit(pfx, mlist(goodi), length(goodi)-2);
  % search for the function minima in the range perturb(2:end)
  %   pf(1)*log10(x).^2+pf(2)*log10(x)+pf(3);
  FF = @(pf, x) polyval(pf, log10(x/p_rng));
  alpha = fminbnd(@(x) FF(pf, x), perturb(min(goodi)), perturb(end));
  % now check how we did
  img.elem_data = x + alpha*dx;
  [dv, opt] = feval(opt.line_search_dv_func, img, data1, N, opt);
  de = feval(opt.line_search_de_func, img, img1, opt);
  meas_err = feval(opt.residual_func, dv, de, W, hp2RtR);
  if opt.verbose > 1
     fprintf('      step size = %0.3g, misfit = %0.3g, expected = %0.3g\n', alpha, meas_err, FF(pf, alpha));
  end

  % check how close we were to the line fit estimate
  % suggest stronger regularization if we're way off
  % (we could automate this update if we wanted to turn on some hueristic behaviour)
  est_err = meas_err / FF(pf, alpha);
  if (opt.verbose > 1) && ((est_err > 1.3) || (est_err < 0.5))
    fprintf('      step misfit missed estimate (%0.1fx est)\n', est_err);
    fprintf('        consider stronger regularization?\n');
  end
else % only two points
  % we provide a FF and pf that will work for plot_line_optimize()
  % this is a straight line between alpha=1 and alpha=1/10
  pf = [];
  FF = @(pf, x)  -(mlist(1) - mlist(end))*0.8*log10(x) + mlist(end);
end
% We save our first choice, in case we are plotting the line search
alpha1 = alpha; % better guess: minima of the fitted curve
meas_err1 = meas_err;

% what if we're making a bad choice?
% if our choice sucked, we've already calculated mlist(:) other choices, go with the minimum
if meas_err > min(mlist)
  [meas_err,mi]= min(mlist);
  alpha= perturb(mi);
  img.elem_data = x + alpha*dx;
  if (length(goodi) > 2) && (opt.verbose > 1)
    fprintf('      did not like our step selection - choose one of perturb values\n');
  end
end

if opt.verbose > 1
   max_alpha_str = '';
   if alpha > perturb(end)-eps
     max_alpha_str = ' (max)';
   end
   fprintf('      step size = %0.3g%s, misfit = %0.3g selected\n', alpha, max_alpha_str, meas_err);
end

% must create plots before changing the perturb values
if opt.line_search_args.plot
  figure;
  plot_line_optimize(perturb, mlist, alpha, meas_err, alpha1, meas_err1, FF, pf);
end

% update perturbations
if meas_err >= mlist(1)
    % this happens when the solution blew up -- the measurement fit was worse than if we did nothing
    if opt.verbose > 1
       fprintf('      reducing perturbations: bad step\n');
    end
    % try a smaller step next time (10x smaller)
    % this keeps the log-space distance between sample points but
    % re-centres around the most recent alpha
    perturb = perturb/10;
else % good step
    if opt.verbose > 1
       fprintf('      update perturbations around step = %0.3g\n', alpha);
    end
    % this keeps the log-space distance between sample points but
    % re-centres around the most recent alpha
    perturb = perturb*(alpha/perturb(end))*2;
end
% jiggle the perturb values by 1% --> if we're stuck in a recursion
% of bad perturb values maybe this is enough to break us out
opt.line_search_args.perturb = perturb .* exp(randn(size(perturb))*0.01);

% Record the corresponding parameters
%img.elem_data= exp(img.logCond);

% we took a bad step, try again but don't recurse indefinitely
if alpha == 0 && retry < 5
  if opt.verbose > 1
     fprintf('    retry#%d (attempt with smaller perturbations)\n', retry+1);
  end
  [alpha, img, dv, opt] = line_search_onm2(imgk, dx, data1, img1, N, W, hp2RtR, dv0, opt, retry+1, pf_max);
end

function plot_line_optimize(perturb, mlist, alpha, meas_err, alpha1, meas_err1, FF, pf)
semilogx(perturb(2:end),mlist(2:end),'xk', 'MarkerSize',10);
hold on;
semilogx(alpha, meas_err,'or', 'MarkerSize',10);
semilogx(alpha1, FF(pf, alpha1), 'ob', 'MarkerSize',10);
semilogx(alpha1, meas_err1, 'xb', 'MarkerSize',10);
% construct the fitted line for plotting
if isnan(perturb(2))
   perturb(2) = perturb(end)/10;
end
p= logspace(log10(perturb(2)),log10(perturb(end)),50);
semilogx(p, FF(pf, p),'k','linewidth',2);
semilogx(p, p*0+mlist(1),'k--','linewidth',1); % alpha=0 should be plotted as well
legend('perturb', 'selected', '1st est', '1st act', 'fit', '\alpha=0');
legend('Location', 'EastOutside');
m = [mlist meas_err meas_err1];
mi=find(isnan(m) | isinf(m)); m(mi) = []; % remove bad values
mr = range(m);
axis([perturb(2) perturb(end) min(m)-mr*0.2 max(m)+mr*0.2]);
xlabel('step size \alpha'); %,'fontsize',20,'fontname','Times')
ylabel('normalized residuals'); %,'fontsize',20,'fontname','Times')
title({sprintf('best alpha = %1.2e',alpha), ...
       sprintf('norm w/o step = %0.4e',mlist(1))}); %,'fontsize',30,'fontname','Times')
%set(gca,'fontsize',20,'fontname','Times');
drawnow; pause(0.5);
