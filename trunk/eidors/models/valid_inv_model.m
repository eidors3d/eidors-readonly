function [pass, err_str] = valid_fwd_model(imdl)
% [pass, err_str] = valid_fwd_model(imdl)
%
% Confirms that a valid forward model structure exists or
% explain why a model is not valid.
%
% If called without an output argument (nargout=0), will
% error out if invalid. Otherwise the function is silent,
% returning an explaination of failures in err_str.

% (C) 2015 Alistair Boyle. License: GPL version 2 or version 3
% $Id$

pass = 1;
err_str = '';

% it's a struct with fields
if ~isstruct(imdl)
   pass = 0;
   err_str = [err_str '- not a struct\n'];
end

% required fields
%      field         type
f = {'name',        'char', ...
     'fwd_model',   'struct', ...
     'solve',       'function', ...
     'reconst_type','char', ... % 'absolute' or 'difference'
     'type',        'char'};
%     'hyperparameter', 'numeric', ... % struct, func, or numeric?
%     'RtR_prior',   'function', ... % function OR numeric
%    'jacobian_bkgnd', 'struct', ... % struct OR numeric
for i=1:length(f)/2
   x=2*(i-1)+1;
   y=x+1;
   if ~isfield(imdl, f{x})
      pass = 0;
      err_str = [err_str '- missing required field: "' f{x} '"\n'];
   elseif strcmp(f{y},'function')
      if ~isfunc(imdl.(f{x}))
         pass = 0;
         err_str = [err_str '- expected function (pointer or string): "' f{x} '"\n'];
      end
   elseif ~isa(imdl.(f{x}), f{y})
      pass = 0;
      err_str = [err_str '- required field "' f{x} '" is not a ' f{y} '\n'];
   end
end
% check the fwd_model
[pass_local, err_str_local] = valid_fwd_model(imdl.fwd_model);
if ~pass_local
   pass = 0;
   disp(err_str_local);
   err_str_local = strrep(err_str_local, ' "', ' "fwd_model.');
   err_str = [err_str err_str_local];
end
clear err_str_local pass_local;
% check for correct 'type'
if ~strcmp(imdl.type, 'inv_model')
   pass = 0;
   err_str = [err_str '- field "type" must be "inv_model"\n'];
end

% optional fields
%      field       type
f = {'rec_model',   'struct'};
for i=1:length(f)/2
   x=2*(i-1)+1;
   y=x+1;
   if isfield(imdl, f{x}) && ~isa(imdl.(f{x}), f{y})
      pass = 0;
      err_str = [err_str '- optional field "' f{x} '" is not a ' f{y} '\n'];
   end
end
% check the rec_model
if isfield(imdl, 'rec_model')
   [pass_local, err_str_local] = valid_fwd_model(imdl.rec_model, 'rec_model');
   if ~pass_local
      pass = 0;
      disp(err_str_local);
      err_str_local = strrep(err_str_local, ' "', ' "rec_model.');
      err_str = [err_str err_str_local];
   end
   clear err_str_local pass_local;
end

% illegal fields
%      field
f = {'inv_model'}; % no recursion
for i=1:length(f)/2
   x=2*(i-1)+1;
   y=x+1;
   if isfield(imdl, f{x})
      pass = 0;
      err_str = [err_str '- illegal field "' f{x} '" found\n'];
   end
end

% result
if ~pass
   err_str = err_str(1:end-2); % drop last \n
end
if ( nargout == 0 ) && ~pass
   error(sprintf(['Reasons:\n' err_str]));
end

function t=isfunc(f)
t=isa(f, 'function_handle') || isa(f, 'char');
