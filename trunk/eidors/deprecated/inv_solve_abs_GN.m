function img= inv_solve_abs_GN(imdl,varargin);
% inv_solve_abs_GN is deprecated in favour of inv_solve_gn

warning('EIDORS:deprecated','INV_SOLVE_ABS_GN is deprecated in favour of INV_SOLVE_GN as of 04-Feb-2016.');
try imdl.inv_solve_gn = imdl.inv_solve_abs_GN; imdl=rmfield(imdl,'inv_solve_abs_GN'); end
imdl = deprecate_imdl_opt(imdl, 'parameters');
imdl = deprecate_imdl_opt(imdl, 'inv_solve');
imdl = deprecate_imdl_opt(imdl, 'inv_solve_core');
img=inv_solve_gn(imdl,varargin{:});
try img.inv_solve_abs_GN = img.inv_solve_gn; img=rmfield(img,'inv_solve_gn'); end

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
      warning('EIDORS:deprecatedParameters',['INV_SOLVE imdl.' opt '.* are deprecated in favour of imdl.inv_solve_gn.* as of 30-Apr-2014.']);
   end

   if ~isfield(imdl, 'inv_solve_gn')
      imdl.inv_solve_gn = imdl.(opt);
   else % we merge
      % merge struct trick from:
      %  http://stackoverflow.com/questions/38645
      for i = fieldnames(imdl.(opt))'
         imdl.inv_solve_gn.(i{1})=imdl.(opt).(i{1});
      end
   end
   imdl = rmfield(imdl, opt);
