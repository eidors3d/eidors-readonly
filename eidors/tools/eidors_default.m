function varargout = eidors_default(varargin)
% EIDORS_DEFAULT - default function handler

optargin = size(varargin,2);

% unit test
if optargin == 1 && ischar(varargin{1}) && strcmp(varargin{1},'UNIT_TEST')
    do_unit_test;
    return;
end

% setters and getters
if optargin > 1 && ischar(varargin{1})
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
            return
    end
end
s = dbstack;
caller  = s(2).name;
varargout{:} = call_default(caller,varargin{:});

function set_default(caller, default)
    global eidors_objects
    eidors_objects.defaults.(caller) = default;


function default = get_default(caller)
    global eidors_objects
    default = eidors_objects.defaults.(caller);

function varargout = call_default(caller,varargin);
default = get_default(caller);
varargout{:} = feval(default,varargin{:});


function do_unit_test
eidors_default('set','eidors_default','test_def');
eidors_default(5);
val = eidors_default(6)
eidors_default('get','eidors_default')