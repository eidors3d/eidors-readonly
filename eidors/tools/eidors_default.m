function varargout = eidors_default(varargin)
% EIDORS_DEFAULT - default function handler

optargin = size(varargin,2);

% unit test
if optargin == 1 && ischar(varargin{1}) && strcmp(varargin{1},'UNIT_TEST')
    do_unit_test;
    return;
end

% setters and getters
if optargin > 0 && ischar(varargin{1})
    switch varargin{1}
        case 'get'
            varargout{1} = get_default(varargin{2});
            return
        case 'set'
            if optargin ~= 3
                error('Wrong number of inputs');
            end
            set_default(varargin{2},varargin{3});
            return
        case 'list'
            varargout{1} = list_defaults;
            return
    end
end
s = dbstack;
caller  = s(2).name;
varargout{:} = call_default(caller,varargin{:});

function list = list_defaults
    global eidors_objects
    try
        list = eidors_objects.defaults;
    catch
        list = [];
    end


function set_default(caller, default)
    global eidors_objects
    eidors_objects.defaults.(caller) = default;


function default = get_default(caller)
    global eidors_objects
    try
        default = eidors_objects.defaults.(caller);
    catch
        error(['No default implementation of ' caller '.']);
    end

function varargout = call_default(caller,varargin);
default = get_default(caller);
varargout{:} = feval(default,varargin{:});


function do_unit_test
fid = fopen('test_def.m','w');
fprintf(fid, 'function out = test_def(in);');
fprintf(fid, 'disp(''test_def'');');
fprintf(fid, 'disp(in);');
fprintf(fid, 'out = in;');
fclose(fid);
eidors_default('set','do_unit_test','test_def');
eidors_default(5);
val = eidors_default(6)
eidors_default('get','do_unit_test')
delete('test_def.m');