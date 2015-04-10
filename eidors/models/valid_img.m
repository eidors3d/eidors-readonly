function pass = valid_img(img)
% [pass, err_str] = valid_img(img)
%
% Confirms that a valid image structure exists or
% explain why an image is not valid.
%
% If called without an output argument (argout=0), will
% error out if invalid. Otherwise the function is silent,
% returning an explaination of failures in err_str.

% (C) 2015 Alistair Boyle. License: GPL version 2 or version 3
% $Id$

pass = 1;
err_str = '';

% it's a struct with fields
if ~isstruct(img)
   pass = 0;
   err_str = [err_str '- not a struct\n'];
end

% required fields
%      field       type
f = {'name',      'char', ...
     'elem_data', 'numeric', ...
     'fwd_model', 'struct', ...
     'type',      'char'};
%     'inv_model', 'struct', ...
for i=1:length(f)/2
   x=2*(i-1)+1;
   y=x+1;
   if ~isfield(img, f{x})
      pass = 0;
      err_str = [err_str '- missing required field: ' f{x} '\n'];
   elseif ~isa(img.(f{x}), f{y})
      pass = 0;
      err_str = [err_str '- required field ' f{x} ' is not a ' f{y}'\n'];
   end
end
% check the fwd_model
[pass_local, err_str_local] = valid_fwd_model(img.fwd_model,'rec_model');
if ~pass_local
   pass = 0;
   disp(err_str_local);
   err_str_local = strrep(err_str_local, ' "', ' "fwd_model.');
   err_str = [err_str err_str_local];
end
clear err_str_local pass_local;
% check for correct 'type'
if ~strcmp(img.type, 'image')
   pass = 0;
   err_str = [err_str '- field "type" must be "image"\n'];
end

% optional fields
%      field       type
f = {'inv_model', 'struct'};
for i=1:length(f)/2
   x=2*(i-1)+1;
   y=x+1;
   if isfield(img, f{x}) && ~isa(img.(f{x}), f{y})
      pass = 0;
      err_str = [err_str '- optional field ' f{x} ' is not a ' f{y}'\n'];
   end
end
% check the inv_model
if isfield(img, 'inv_model')
   [pass_local, err_str_local] = valid_inv_model(img.inv_model);
   if ~pass_local
      pass = 0;
      disp(err_str_local);
      err_str_local = strrep(err_str_local, ' "', ' "inv_model.');
      err_str = [err_str err_str_local];
   end
   clear err_str_local pass_local;
end

% illegal fields
%      field       type
f = {'imdl', 'inv_mdl', 'fmdl', 'fwd_mdl'};
for i=1:length(f)
   x=i;
   if isfield(img, f{x})
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
