function [img fmin res] = line_optimize(imgk, dx, data0, opt)
%LINE_OPTIMIZE Cheap line optimizer
% [img fmin res] = line_optimize(imgk, dx, data0, opt)
% img     : output image
% fmin    : optimal step size
% res     : value of the objective function
% imgk    : starting image
% dx      : step direction
% data0   : data to fit
% opt     : options structure
%
% Options:
%   opt.perturb  : vector of step sizes to try (default: [0.1 0.5 1.0])
%   opt.min_value: lower limit of img values (default: -Inf)
%   opt.max_value: upper limit of img values (default: Inf)
%   opt.objective_function: @my_objective_fun
%       handle to an objective funtion with the following signature
%           val = my_objective_fun(data0, data, img0, img, opt)
%       where:
%           data0   : data to fit
%           data    : data simulated on the current image
%           img0    : starting image
%           img     : current image
%           opt     : the entire option structure (can be used to pass
%           additional parameters)
%       The default objective function is:
%           val = norm(calc_difference_data(data0,data, img0.fwd_model));
%
% Note that the value of fmin will be limited by the maximum and minimum
% given in the opt.perturb vector.
%
% See also: INV_SOLVE_ABS_GN

% (C) 2010-2013 Copyright Bartlomiej Grychtol, Andy Adler & Nolwenn Lesparre.
% License: GPL version 2 or 3.
% $Id$


if nargin <4 || isempty(opt)
   opt = struct;
end

opt = parse_options(opt);

img = imgk;
for i = 1:length(opt.perturb)
   img = calc_perturb(imgk,opt.perturb(i),dx, opt);
   vsim = fwd_solve(img);
   mlist(i) = feval(opt.objective_func,data0,vsim,imgk,img,opt);
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
img  = calc_perturb(imgk,fmin,dx,opt);
vsim = fwd_solve(img);
res  = feval(opt.objective_func,data0,vsim,imgk,img,opt);



function img = calc_perturb(imgk, p, dx, opt)
   imgk = physics_data_mapper(imgk);
   img = imgk;
   img.elem_data = imgk.elem_data + p*dx;
   img = apply_limits(img,opt);
   img = physics_data_mapper(img, 1);


function img = apply_limits(img,opt)
   img.elem_data(img.elem_data > opt.max_value) = opt.max_value;
   img.elem_data(img.elem_data < opt.min_value) = opt.min_value;


function val = default_obj_fun(data0, data, img0, img, opt)
dv = calc_difference_data(data0,data,img0.fwd_model);
val = norm(dv);


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

if ~isfield(opt,'objective_func')
    opt.objective_func = @default_obj_fun;
end
