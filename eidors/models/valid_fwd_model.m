function [pass, err_str] = valid_fwd_model(fmdl, type)
% [pass, err_str] = valid_fwd_model(fmdl, [type])
%
% Confirms that a valid forward model structure exists or
% explain why a model is not valid.
%
% The model is assumed to be a fwd_model and all fields are
% checked unless type='rec_model'. A reconstruction model
% (rec_model) is only checked to confirm it has a valid mesh
% associated with it.
%
% If called without an output argument (nargout=0), will
% error out if invalid. Otherwise the function is silent,
% returning an explaination of failures in err_str.

% (C) 2015 Alistair Boyle. License: GPL version 2 or version 3
% $Id$

pass = 1;
err_str = '';

if nargin < 2
   type = 'fwd_model';
end
if ~any(strcmp(type, {'fwd_model','rec_model'}))
   error('unexpected "type"');
end

% it's a struct with fields
if ~isstruct(fmdl)
   pass = 0;
   err_str = [err_str '- not a struct\n'];
end

% required fields
%      field         type
f = {'name',        'char', ...
     'nodes',       'numeric', ... % uintX
     'elems',       'numeric', ... % uintX
     'gnd_node',    'numeric', ... % uintX
     'electrode',   'struct', ...
     'stimulation', 'struct', ...
     'solve',       'function', ...
     'system_mat',  'function', ...
     'jacobian',    'function', ...
     'normalize_measurements', 'numeric', ... % logical
     'type',        'char'};
%     'boundary',   'numeric', ... % uintX OPTIONAL
%     'meas_select','numeric', ... % uintX OPTIONAL
% reduced set of requirements for a rec_model
if strcmp(type, 'rec_model')
   f = {'name',     'char', ...
        'nodes',    'numeric', ...
        'elems',    'numeric'};
end
for i=1:length(f)/2
   x=2*(i-1)+1;
   y=x+1;
   if ~isfield(fmdl, f{x})
      pass = 0;
      err_str = [err_str '- missing required field: "' f{x} '"\n'];
   elseif strcmp(f{y},'function')
      if ~isfunc(fmdl.(f{x}))
         pass = 0;
         err_str = [err_str '- expected function (pointer or string): "' f{x} '"\n'];
      end
   elseif ~isa(fmdl.(f{x}), f{y})
      pass = 0;
      err_str = [err_str '- required field "' f{x} '" is not a ' f{y} '\n'];
   end
end
if ~strcmp(type, 'rec_model')
   % validate electrode struct
   nn = size(fmdl.nodes,1); % number of nodes
   for i=1:length(fmdl.electrode)
      [pass, err_str] = valid_elec(fmdl.electrode(i), i, nn, pass, err_str);
   end
   % validate stimulation struct
   nel = length(fmdl.electrode); % number of electrodes
   for i=1:length(fmdl.stimulation)
      [pass, err_str] = valid_stim(fmdl.stimulation(i), i, nel, pass, err_str);
   end
end
% check for correct 'type'
if ~any(strcmp(fmdl.type, {'fwd_model', 'rec_model'}))
   pass = 0;
   err_str = [err_str '- field "type" must be "fwd_model" or "rec_model"\n'];
end

% optional fields
%      field       type
f = {'boundary',   'numeric', ... % uintX
     'meas_select','logical'};
for i=1:length(f)/2
   x=2*(i-1)+1;
   y=x+1;
   if isfield(fmdl, f{x}) && ~isa(fmdl.(f{x}), f{y})
      pass = 0;
      err_str = [err_str '- optional field "' f{x} '" is not a ' f{y} '\n'];
   end
end

% illegal fields (common typos, etc)
%      field
f = {'fwd_model', ...
     'rec_model', ... % no recursion
     'stim', ... % not short form
     'elec', ...
     'fmdl', ...
     'fwd_mdl', ...
     'cmdl', ...
     'rmdl', ...
     'rec_mdl', ...
     'electrodes', ... % not plural
     'stimulations', ...
     'node', ... % plural
     'elem'};
for i=1:length(f)
   x=i;
   if isfield(fmdl, f{x})
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

function [pass, err_str] = valid_stim(stim, i, nel, pass, err_str)
pass_local = 1;
% required fields
%      field         type
f = {'stim_pattern', 'numeric', ...
     'meas_pattern', 'numeric'};
for i=1:length(f)/2
   x=2*(i-1)+1;
   y=x+1;
   if ~isfield(stim, f{x})
      pass_local = 0;
      err_str = [err_str '- missing required field: "stimulation(' num2str(i) ').' f{x} '"\n'];
   elseif ~isa(stim.(f{x}), f{y})
      pass_local = 0;
      err_str = [err_str '- required field "stimulation(' num2str(i) ').' f{x} '" is not a ' f{y} '\n'];
   end
end
% optional fields
%      field       type
f = {'stimulation',  'char'};
for i=1:length(f)/2
   x=2*(i-1)+1;
   y=x+1;
   if isfield(stim, f{x}) && ~isa(stim.(f{x}), f{y})
      pass_local = 0;
      err_str = [err_str '- optional field "stimulation(' num2str(i) ').' f{x} '" is not a ' f{y} '\n'];
   end
end
% we don't need to check further
if ~pass_local
   pass = 0;
   return
end
clear pass_local;
% we expect only stimulation, meas/stim_pattern fields
len = length(fields(stim));
if (len > 3) || ((len > 2) && ~any(strcmp(fields(stim),'stimulation')))
   pass = 0;
   err_str = [err_str '- extraneous fields in "stimulation(' num2str(i) ')\n'];
end
if size(stim.stim_pattern,1) > nel
   pass = 0;
   err_str = [err_str '- field "stimulation(' num2str(i) ').stim_pattern has more rows than electrodes in the model\n'];
end
if size(stim.meas_pattern,2) > nel
   pass = 0;
   err_str = [err_str '- field "stimulation(' num2str(i) ').meas_pattern has more columns than electrodes in the model\n'];
end

function [pass, err_str] = valid_elec(elec, i, nn, pass, err_str);
pass_local = 1;
% required fields
%      field         type
f = {'nodes',     'numeric', ...
     'z_contact', 'numeric'};
for i=1:length(f)/2
   x=2*(i-1)+1;
   y=x+1;
   if ~isfield(elec, f{x})
      pass_local = 0;
      err_str = [err_str '- missing required field: "electrode(' num2str(i) ').' f{x} '"\n'];
   elseif ~isa(elec.(f{x}), f{y})
      pass_local = 0;
      err_str = [err_str '- required field "electrode(' num2str(i) ').' f{x} '" is not a ' f{y} '\n'];
   end
end
% we don't need to check further
if ~pass_local
   pass = 0;
   return
end
clear pass_local;
% check that 'nodes' are valid
if any((elec.nodes > nn) | (elec.nodes < 1))
   pass = 0;
   err_str = [err_str '- field "electrode(' num2str(i) ').nodes do not exist on the model.nodes\n'];
end
