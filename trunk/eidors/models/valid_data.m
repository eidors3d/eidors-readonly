function pass = valid_data(data)
% [pass, err_str] = valid_data(data)
%
% Confirms that a valid measurement data structure exists or
% explain why an structure is not valid.
%
% If called without an output argument (argout=0), will
% error out if invalid. Otherwise the function is silent,
% returning an explaination of failures in err_str.

% (C) 2015 Alistair Boyle. License: GPL version 2 or version 3
% $Id$

pass = 1;
err_str = '';

% it *can* be just straight up numbers
if isnumeric(data)
   return;
end

% it's a struct with fields
if ~isstruct(data)
   pass = 0;
   err_str = [err_str '- not a struct or numeric data\n'];
end

% required fields
%      field       type
f = {'name',      'char', ...
     'meas',      'numeric', ...
     'time',      'numeric', ... % OPTIONAL?
     'type',      'char'};
for i=1:length(f)/2
   x=2*(i-1)+1;
   y=x+1;
   if ~isfield(data, f{x})
      pass = 0;
      err_str = [err_str '- missing required field: ' f{x} '\n'];
   elseif ~isa(data.(f{x}), f{y})
      pass = 0;
      err_str = [err_str '- required field ' f{x} ' is not a ' f{y}'\n'];
   end
end
% check for correct 'type'
if ~strcmp(data.type, 'data')
   pass = 0;
   err_str = [err_str '- field "type" must be "data"\n'];
end

% optional fields
%      field       type
f = {};
for i=1:length(f)/2
   x=2*(i-1)+1;
   y=x+1;
   if isfield(data, f{x}) && ~isa(data.(f{x}), f{y})
      pass = 0;
      err_str = [err_str '- optional field ' f{x} ' is not a ' f{y}'\n'];
   end
end

% illegal fields
%      field
f = {'maes'};
for i=1:length(f)
   x=i;
   if isfield(data, f{x})
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
