function [img fmin res] = line_optimize(imgk, dx, data1, opt)

if nargin <4 || isempty(opt)
   opt = struct;
end

opt = parse_options(opt);

imgk = physics_data_mapper(imgk);

img = imgk;
for i = 1:length(opt.perturb)
   img = calc_perturb(imgk,opt.perturb(i),dx, opt);
   vsim = fwd_solve(img);
   % aren't we better off just calculating the difference ourselves?
   dv = calc_difference_data( vsim, data1, img.fwd_model);
   mlist(i) = norm(dv);
end

% perform quadratic line fit in log space
p10 = log10(opt.perturb);

pf = polyfit(p10, mlist, 2);
fmin = -pf(2)/pf(1)/2; % poly minimum for a 2nd order poly
val = polyval(pf, fmin);

if val > min(mlist) % fit didn't find the minimum, correct
   [mlist_o ik] = sort(mlist);
   % in case the user provided a 0 perturbation, avoid it
   if opt.perturb( ik(1)) ~= 0
      fmin = opt.perturb(1);
   else
      fmin = opt.perturb(2);
   end
end

fmin = 10^fmin; % convert back to linear

% limit to the values given in perturb
if fmin < min( opt.perturb )
   fmin = min( opt.perturb );
end

if fmin > max( opt.perturb )
   fmin = max( opt.perturb );
end

% RETURN VALUES
img  = calc_perturb(imgk,opt.perturb(i),dx,opt);
vsim = fwd_solve(img);
res  = calc_difference_data( vsim, data1, img.fwd_model);



function img = calc_perturb(imgk, p, dx, opt)
   img = imgk;
   img.elem_data = imgk.elem_data + p*dx;
   img = apply_limits(img,opt);
   img = physics_data_mapper(img, 1);


function img = apply_limits(img,opt)
img.elem_data(img.elem_data > opt.max_value) = opt.max_value;
img.elem_data(img.elem_data < opt.min_value) = opt.min_value;


function opt = parse_options(opt)
if ~isfield(opt,'perturb')
   opt.perturb = [0.1 0.5 1.0];
end

if ~isfield(opt,'min_value')
   opt.min_value = -Inf;
end


if ~isfield(opt,'max_value')
   opt.max_value = Inf;
end