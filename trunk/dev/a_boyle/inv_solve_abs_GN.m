function img= inv_solve_abs_GN( inv_model, data1);
%function img= inv_solve_abs_GN( inv_model, data1);
% INV_SOLVE_ABS_GN
% This function calls INV_SOLVE_ABS_CORE to find a Gauss-Newton
% iterative solution.
%
% img = inv_solve_abs_GN( inv_model, data1 )
%   img        => output image data (or vector of images)
%   inv_model  => inverse model struct
%   data1      => measurements
%
% See INV_SOLVE_ABS_CORE for arguments, options and parameters.
%
% (C) 2014 Alistair Boyle
% License: GPL version 2 or version 3
% $Id$

%--------------------------
% UNIT_TEST?
if isstr(inv_model) && strcmp(inv_model,'UNIT_TEST'); img = do_unit_test; return; end

%if isfield(inv_model, 'inv_solve_abs_core') && isfield(inv_model, 'inv_solve_abs_GN')
%   error('merging inv_solve_abs_core and inv_solve_abs_GN options is broken, please use only one method');
%end
%% merge legacy options locations
inv_model = deprecate_imdl_opt(inv_model, 'parameters');
inv_model = deprecate_imdl_opt(inv_model, 'inv_solve');
inv_model = deprecate_imdl_opt(inv_model, 'inv_solve_abs_core');

% inv_model.inv_solve_abs_GN -> inv_solve_abs_core
if isfield(inv_model, 'inv_solve_abs_GN')
   if isfield(inv_model, 'inv_solve_abs_core')
      error('inv_model.inv_solve_abs_GN replaces inv_model.inv_solve_abs_core, parameters will be lost');
   end
   inv_model.inv_solve_abs_core = inv_model.inv_solve_abs_GN;
   inv_model = rmfield(inv_model, 'inv_solve_abs_GN');
end

img = inv_solve_abs_core(inv_model, data1);

if isfield(img, 'inv_solve_abs_core')
  img.inv_solve_abs_GN = img.inv_solve_abs_core;
  img=rmfield(img, 'inv_solve_abs_core');
end

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
      warning('EIDORS:deprecatedParameters',['INV_SOLVE inv_model.' opt '.* are deprecated in favor of inv_model.inv_solve_abs_GN.* as of 30-Apr-2014.']);
   end

   if ~isfield(imdl, 'inv_solve_abs_GN')
      imdl.inv_solve_abs_GN = imdl.(opt);
   else % we merge
      % merge struct trick from:
      %  http://stackoverflow.com/questions/38645
      for i = fieldnames(imdl.(opt))'
         imdl.inv_solve_abs_GN.(i{1})=imdl.(opt).(i{1});
      end
   end
   imdl = rmfield(imdl, opt);

function pass = do_unit_test()
   pass = inv_solve_abs_core('UNIT_TEST', 'inv_solve_abs_GN');
