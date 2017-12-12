function  [alpha, img, dv, opt] = line_search_onm2(imgk, dx, data1, img1, N, W, hps2RtR, hpt2LLt, dv0, opt, retry, pf_max)
% function  [alpha, img, dv, opt] = line_search_onm2(imgk, dx, data1, img1, N, W, hps2RtR, hpt2LLt, dv0, opt)
% line search function with a fitted polynomial of O(n-2) where n is the number of perturbations
% (C) 2013 Alistair Boyle
% License: GPL version 2 or version 3

perturb = calc_perturb(imgk, dx, opt);

if nargin < 11
  retry = 0;
end

if nargin < 12
  pf_max = length(perturb)-2;
end

% fwd_solve is the most expensive part generally, count how many we do
if ~isfield(opt, 'fwd_solutions')
   opt.fwd_solutions = 0;
end
x = imgk.elem_data;

if(perturb(1) ~= 0)
  error('first perturbation min(inv_model.inv_solve_abs_{GN,CG,core}.line_search.perturb) expects alpha=0');
end

% Compute the forward model for each perturbation step
img = imgk;
% mlist is our search result for each alpha value, perturb(i)
%  -- NaN: initiailized but not calculated
%  -- -Inf: should not occur we have code that converts calculated NaNs and -Inf to +Inf
%  -- +Inf: calculated value was bad, ignore it
if opt.verbose > 1
   fprintf('      i     = ');
   fprintf('    [%d]  \t', 1:length(perturb));
   fprintf('\n');
   fprintf('      alpha = ');
   fprintf(' %8.3g\t', perturb);
   fprintf('\n');
   fprintf('              ');
end
mlist= ones(size(perturb))*NaN; % init
for i = 1:length(perturb)
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
      mlist(i) = feval(opt.residual_func, dv, de, W, hps2RtR, hpt2LLt);
      if any(isnan(mlist(i)) | isinf(mlist(i)))
         mlist(i) = +Inf; % NaN or Inf are converted to Inf, since we use NaN to indicate initialized but not calculated
      end
   end
   if opt.verbose > 1
      fprintf(' %8.3g\t',mlist(i));
   end
   if mlist(i)/mlist(1) > 1e10
      if opt.verbose > 1
         for j=(i+1):length(perturb)
            fprintf('   [skip]\t');
         end
      end
      break
   end
end
if opt.verbose > 1
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
  meas_err = feval(opt.residual_func, dv, de, W, hps2RtR, hpt2LLt);
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
  dv = []; % we'll need to recalculate this later since we didn't keep it
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
% 
% keyboard
% must create plots before changing the perturb values
if opt.line_search_args.plot
   clf;
   plot_line_optimize(perturb, mlist, alpha, meas_err, alpha1, meas_err1, FF, pf);
   if isfield(opt,'fig_prefix') % TODO assign from base options if set
      k=1; % TODO iteration count; TODO retry count
      print('-dpdf',sprintf('%s-ls%d-retry%d',opt.fig_prefix,k,retry));
      print('-dpng',sprintf('%s-ls%d-retry%d',opt.fig_prefix,k,retry));
      saveas(gcf,sprintf('%s-ls%d-retry%d.fig',opt.fig_prefix,k,retry));
   end
end

% update perturbations
if meas_err >= mlist(1)
    if mlist(1)*1.05 < mlist(goodi(end))
       % this happens when the solution blew up -- the measurement fit was worse than if we did nothing
       if opt.verbose > 1
          fprintf('      reducing perturbations /10: bad step\n');
       end
       % try a smaller step next time (10x smaller)
       % this keeps the log-space distance between sample points
       perturb = perturb/10;
    elseif perturb(end) > 1.0-10*eps
       if opt.verbose > 1
          fprintf('      expanding perturbations x10: ... but we''d be searching past alpha=1.0, giving up\n');
       end
       return % we give up early
    elseif perturb(end)*10 > 1.0-10*eps
       if opt.verbose > 1
          fprintf('      expanding perturbations (limit alpha=1.0): bad step\n');
       end
       perturb = perturb/perturb(end); % ... max(alpha)=1.0
    else % we didn't really get any difference in solutions
       % this happens when the perturbations are too small, we are too close to
       % the current solution
       if opt.verbose > 1
          fprintf('      expanding perturbations x10: bad step\n');
       end
       % try a larger step next time (10x larger)
       % this keeps the log-space distance between sample points
       perturb = perturb*10;
    end
else % good step
    % stretch out the perturbations if we're not making much progress
    if all(mlist(goodi)/mlist(1)-1 > -10*opt.dtol) && ...
       (perturb(end)*10 < 1.0+10*eps)
       if opt.verbose > 1
          fprintf('      expand perturbations x10 for next iteration\n');
          fprintf('      (didn''t make much progress this iteration)\n');
       end
       opt.line_search_args.perturb = opt.line_search_args.perturb*10;
    else % or just recentre around our best answer
       % this keeps the log-space distance between sample points but
       % re-centres around the most recent alpha
       if opt.verbose > 1
          fprintf('      update perturbations around step = %0.3g (limit alpha=1.0)\n', alpha);
       end
       if alpha/perturb(end)*2 > 1.0 - 10*eps
          perturb = perturb/perturb(end);
       else
          perturb = perturb*(alpha/perturb(end))*2;
       end
    end
end
% jiggle the perturb values by 1% --> if we're stuck in a recursion
% of bad perturb values maybe this is enough to break us out
perturb = perturb .* exp(randn(size(perturb))*0.01);
% fix if we exceeded alpha=1.0
if perturb(end) > 1.0 - eps
   perturb = perturb/perturb(end);
end
opt.line_search_args.perturb = perturb;

% Record the corresponding parameters
%img.elem_data= exp(img.logCond);

% we took a bad step, try again but don't recurse indefinitely
if alpha == 0 && retry < 5
  if opt.verbose > 1
     fprintf('    retry#%d (attempt with new perturbations)\n', retry+1);
  end
  [alpha, img, dv, opt] = line_search_onm2(imgk, dx, data1, img1, N, W, hps2RtR, hpt2LLt, dv0, opt, retry+1, pf_max);
end

%%%
% calculate search values for \alpha
% 1. sort in ascending order
% 2. scale to a range within finite precision (don't waste time searching when
%    results will be inf/nan)
function perturb = calc_perturb(imgk, dx_in, opt)
if opt.verbose > 1
   disp('      line search (finite precision) limits');
end
% scale
% When log numbers are converted to base numbers, they frequently result in Inf
% when dx is large.  We scale 'perturb' here, so that we are line searching
% within a numerically stable region for our finite precision numbers.
%   log(realmax) = 709.7827
% log10(realmax) = 308.2547
% - canonicalize the img data, so we don't have to deal with default forms
if ~isfield(imgk, 'current_params')
   imgk.current_params = {'conductivity'};
end
if ~isfield(imgk, 'params_sel')
   imgk.params_sel = {[1:length(imgk.elem_data)]};
end
if ~iscell(imgk.current_params)
   imgk.current_params = {imgk.current_params};
end
if ~iscell(imgk.params_sel)
   imgk.params_sel = {imgk.params_sel};
end

if isfield(imgk, 'inv_model') && isfield(imgk.inv_model, 'fwd_model')
   md = max(range(imgk.inv_model.fwd_model.nodes)); % model range (coordinates)
end

% determine the maximum alpha
max_alpha = +inf;
min_alpha = +inf;
for i=1:length(imgk.current_params)
   p = imgk.current_params{i};
   s = imgk.params_sel{i};
   x = imgk.elem_data(s);
   dx = dx_in(s);
   % TODO these could be based on the limits provided as args to inv_solve_abs_GN, instead they are hardcoded here...
   is_mvmt = (length(p) >= 8) && strcmp(p(end-7:end),'movement');
   if strcmp(p(1:4), 'log_')
      lp = log(realmax/2); % largest positive floating point number (double): Limit_Positive
      ln = -inf; % or = log(realmin/2); % largest negative floating point number (double): Limit Negative
      % for log space, we should have an ln = -inf --> exp(-900) = 0
      if is_mvmt
         lp = log(md);
      end
   elseif strcmp(p(1:6), 'log10_')
      lp = log10(realmax/2); % largest positive floating point number (double): Limit_Positive
      ln = -inf; % or = log10(realmin/2); % largest negative floating point number (double): Limit Negative
      % for log10 space, we should have an ln = -inf --> 10.^-900 = 0
      if is_mvmt
         lp = log10(md);
      end
   else
      lp = +realmax/2;
      ln = -realmax/2;
      if is_mvmt
         lp = +md;
         lp = -md;
      end
   end
   % lower limit on \alpha prior to x = x + alpha*dx --> +/-inf; % (explodes)
   %   \alpha_min = ((max or min) - x) / \delta_x
   au=(lp-x)./dx; au(dx<=0)=NaN; au(isnan(au))=+inf; au=min(au);
   a_max = au;
   au=(ln-x)./dx; au(dx>=0)=NaN; au(isnan(au))=+inf; au=min(au);
   if (au < a_max)
      a_max = au;
   end
   % lower limit on \alpha prior to x == x + alpha*dx; % (no change)
   %   \alpha_min = \epsilon / \delta_x
   if is_mvmt
      al=1e-3./abs(dx); % don't care about movement less than 1mm
      % TODO configurable? 'reconstruction tolerance'?
   else
      al=eps(x)./abs(dx);
   end
   al(isinf(al))=NaN; al(isnan(al))=realmax; al=min(al);
   if isnan(al)
      a_min = 0;
   else
      a_min = al;
   end
   if opt.verbose > 1
      fprintf('        %s: alpha range = %0.3g -- %0.3g\n', p, a_min, a_max);
   end
   % adjust global limits
   if a_max < max_alpha
      max_alpha = a_max;
   end
   if a_min < min_alpha
      min_alpha = a_min;
   end
end

% sort
p=sort(opt.line_search_args.perturb);
% scale
if (p(end) > max_alpha) || (p(2) < min_alpha)
   p(p<realmin/2) = [];
   p=log10(p); ap=log10(max_alpha); an=log10(min_alpha);
   if range(p) >  ap-an
      p=p*(ap-an)/range(p);
   end
   if p(end) > ap
      p=p-(max(p)-ap);
   elseif p(1) < an
      p=p+(an-min(p));
   end
   p=[0 10.^p];
   if opt.verbose > 1
      fprintf('        alpha (before) = ');
      fprintf('%0.3g ', sort(opt.line_search_args.perturb));
      fprintf('\n');
      fprintf('        alpha (after)  = ');
      fprintf('%0.3g ', p);
      fprintf('\n');
   end
else
   if opt.verbose > 1
      fprintf('        alpha (unchanged) = ');
      fprintf('%0.3g ', p);
      fprintf('\n');
   end
end
perturb=p;


%%% plot the line optimization results
% 1. search locations
% 2. line fit
% 3. selected minima and test point result
% 4. selected \alpha
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
if mr < max(m)*1e-14
   mr = 1e-14;
end
axis([perturb(2) perturb(end) min(m)-mr*0.2 max(m)+mr*0.2]);
xlabel('step size \alpha'); %,'fontsize',20,'fontname','Times')
ylabel('normalized residuals'); %,'fontsize',20,'fontname','Times')
title({sprintf('best alpha = %1.2e',alpha), ...
       sprintf('norm w/o step = %0.4e',mlist(1))}); %,'fontsize',30,'fontname','Times')
%set(gca,'fontsize',20,'fontname','Times');
drawnow;

function x=range(y)
x=max(y)-min(y);
