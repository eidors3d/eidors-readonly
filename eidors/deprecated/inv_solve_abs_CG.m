function img= inv_solve_abs_CG(imdl,varargin);
% inv_solve_abs_CG is deprecated in favour of inv_solve_cg

warning('EIDORS:deprecated','INV_SOLVE_ABS_CG is deprecated in favour of INV_SOLVE_CG as of 04-Feb-2016.');
try imdl.inv_solve_cg = imdl.inv_solve_abs_CG; imdl=rmfield(imdl,'inv_solve_abs_CG'); end
% merge legacy options locations
imdl = deprecate_imdl_opt(imdl, 'parameters');
imdl = deprecate_imdl_opt(imdl, 'inv_solve');
imdl = deprecate_imdl_opt(imdl, 'inv_solve_core');
img=inv_solve_cg(imdl,varargin{:});
try img.inv_solve_abs_CG = img.inv_solve_cg; img=rmfield(img,'inv_solve_cg'); end

function imdl = deprecate_imdl_opt(imdl,opt)
   if ~isfield(imdl, opt)
      return;
   end
   if ~isstruct(imdl.(opt))
      error(['unexpected imdl.' opt ' where ' opt ' is not a struct... i do not know what to do']);
   end

   % warn on anything but imdl.inv_solve.calc_solution_error
   Af = fieldnames(imdl.(opt));
   if ~strcmp(opt, 'inv_solve') || (length(Af(:)) ~= 1) || ~strcmp(Af(:),'calc_solution_error')
      disp(imdl)
      disp(imdl.(opt))
      warning('EIDORS:deprecatedParameters',['INV_SOLVE imdl.' opt '.* are deprecated in favor of imdl.inv_solve_cg.* as of 30-Apr-2014.']);
   end

   if ~isfield(imdl, 'inv_solve_cg')
      imdl.inv_solve_cg = imdl.(opt);
   else % we merge
      % merge struct trick from:
      %  http://stackoverflow.com/questions/38645
      for i = fieldnames(imdl.(opt))'
         imdl.inv_solve_cg.(i{1})=imdl.(opt).(i{1});
      end
   end
   imdl = rmfield(imdl, opt);
