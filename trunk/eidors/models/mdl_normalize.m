function out = mdl_normalize(mdl, val)
%MDL_NORMALIZE Check or set the normalize_measurements flag on a model.
%  OUT = MDL_NORMALIZE(MDL) returns the value of the 
%  'normalize_measurements' or 'normalize' field on MDL. If absent, the 
%  default value is returned as set with EIDORS_DEFAULT. 
%
%  MDL = MDL_NORMALIZE(MDL, VAL) sets the 'normalize_measurements' to VAL
%
% See also: FWD_MODEL_PARAMETERS, EIDORS_DEFAULT, EIDORS_STARTUP

% (C) 2012 Bartlomiej Grychtol. License: GPL v2 or v3
% $Id$

if ischar(mdl) && strcmp(mdl, 'UNIT_TEST'); do_unit_test; return; end

% check that we got a fwd_model
try 
    if ~strcmp(mdl.type,'fwd_model'), error;end;
catch
    error('The normalize_measurements flag must be on a fwd_model');
end
% do the real thing
switch nargin
    case 1
        out = get_flag(mdl);
        if isempty(out)
            out = feval('eidors_default');
            warning('EIDORS:normalize_flag_not_set', ...
                'normalize_measurements flag not set on model. Using default: %d',out);
        end
    case 2
        out = set_flag(mdl,val);
    otherwise
        error('Wrong number of parameters');
end

function out = get_flag(mdl)
 out = [];
% iterate over fields
ff = fieldnames(mdl);
for i = 1:length(ff);
    % if both are set, you're in trouble
    if any(strcmp(ff{i},{'normalize','normalize_measurements'}))
        out = mdl.(ff{i});
    else
        % won't hit 'normal' or anything that starts with 'normals'. 
        token = regexp(ff{i},'normal[^s]+.*','match','once'); 
        if ~isempty(token);
            warning('Suspect field "%s" found on the model',token);
        end
        
    end
end

function mdl = set_flag(mdl, val) 
mdl.normalize_measurements = val;

function do_unit_test
s = warning('query', 'backtrace');
warning off backtrace;
eidors_default('set','mdl_normalize',@(x) 0);
mdl = mk_common_model('c2t2',16);
fmdl = mdl.fwd_model;
eidors_msg('Does it come normalized?');
eidors_msg(mdl_normalize(fmdl));
eidors_msg('Normalize it');
fmdl = mdl_normalize(fmdl,1);
eidors_msg('How about now?');
eidors_msg(mdl_normalize(fmdl));
eidors_msg('De-normalize it');
fmdl = mdl_normalize(fmdl,0);
eidors_msg('And now?');
eidors_msg(mdl_normalize(fmdl));
eidors_msg('Remove the flag');
fmdl = rmfield(fmdl,'normalize_measurements');
eidors_msg(mdl_normalize(fmdl));
eidors_msg('And now let''s set it to something stupid');
fmdl.normalise_measurements = 1; % s, throw a warning, return default
eidors_msg(mdl_normalize(fmdl));
fmdl.normalize_my_gran = 1; % throw a warning, return default
eidors_msg(mdl_normalize(fmdl));
fmdl.normalize = 1; % accept happily
eidors_msg(mdl_normalize(fmdl));
eidors_msg('But an imdl is not acceptable!');
try
eidors_msg(mdl_normalize(mdl));
catch
    eidors_msg('Just so.');
end
warning(s);
