function varargout = eidors_default(varargin)
%EIDORS_DEFAULT Default function handler.
%  EIDORS_DEFAULT redirects calls to the default implementation of the 
%  calling function's functionality, e.g. if called from CALC_JACOBIAN it 
%  will call JACOBIAN_ADJOINT. The default implementations are set at 
%  startup in EIDORS_STARTUP, and can be manipulated as follows:
%
%  eidors_default('set','caller','implementation')  - set the default
%  implementation of 'caller' function to 'implementation'.
%
%  eidors_default('get','caller') - get the current default implementation
%  for 'caller'
% 
%  eidors_default('list') - list all stored defaults
%
%  NOTE: This is an internal function not meant to be called directly other
%  than as indicated above. Rather, where the algorithm to be used by a 
%  certain function is passed to it as a field on an eidors model
%  structure, that field can be set to 'eidors_default'. For example, 
%  the inverse solver to be used is in inv_solve is specified in 
%  imdl.solve. If imdl.solve = 'eidors_default', a call to inv_solve(imdl)
%  will use the default implementation.

% (C) Bartlomiej Grychtol, 2012. License: GPL v. 2 or 3.
% $Id$

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
varargout{1:nargout} = call_default(caller,varargin{:});

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
varargout{1:nargout} = feval(default,varargin{:});


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
