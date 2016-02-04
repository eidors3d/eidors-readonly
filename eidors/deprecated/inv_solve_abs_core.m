function img= inv_solve_abs_core(imdl,varargin);
% inv_solve_abs_core is deprecated in favour of inv_solve_core

warning('EIDORS:deprecated','INV_SOLVE_ABS_CORE is deprecated in favour of INV_SOLVE_CORE as of 04-Feb-2016.');
try imdl.inv_solve_core = imdl.inv_solve_abs_core; imdl=rmfield(imdl,'inv_solve_abs_core'); end
imdl = deprecate_imdl_opt(imdl, 'parameters');
imdl = deprecate_imdl_opt(imdl, 'inv_solve');
imdl = deprecate_imdl_opt(imdl, 'inv_solve_core');
img=inv_solve_core(imdl,varargin{:});
try img.inv_solve_abs_core = img.inv_solve_core; img=rmfield(img,'inv_solve_core'); end

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
      warning('EIDORS:deprecatedParameters',['INV_SOLVE inv_model.' opt '.* are deprecated in favor of inv_model.inv_solve_core.* as of 30-Apr-2014.']);
   end

   if ~isfield(imdl, 'inv_solve_core')
      imdl.inv_solve_core = imdl.(opt);
   else % we merge
      % merge struct trick from:
      %  http://stackoverflow.com/questions/38645
      for i = fieldnames(imdl.(opt))'
         imdl.inv_solve_core.(i{1})=imdl.(opt).(i{1});
      end
   end
   imdl = rmfield(imdl, opt);
