function out = eidors_debug(command, fstr)
%EIDORS_DEBUG Global managment of debug flags
% eidors_debug('on','function') 
%     switches ON the debug flag for string constant 'function'
% eidors_debug('off','function') 
%     switches OFF the debug flag for string constant 'function'
% eidors_debug('off','all') 
%     switches OFF the debug flag for all functions
% eidors_debug('query','function') 
%     returns the current state of the debug flag for 'function'
% 
% Note that 'function' can be any string constant (other than 'all'), but 
% it is recommended to be an actual eidors function name. To enable 
% debugging on selected parts of the code, follow this convention:
%    eidors_function:part1
%    eidors_function:subfunction
%    eidors_function:subfunction:part1
%
% Example use in code (my_function.m):
% if eidors_debug('query','my_function')
%     disp(val)
% end

% (C) 2013 Bartlomiej Grychtol.
% License: GPL version 2 or 3
% $Id$



if ischar(command) && strcmp(command, 'UNIT_TEST'), do_unit_test; end

global eidors_objects;

if ~isfield(eidors_objects, 'debug_enabled_on')
   eidors_objects.debug_enabled_on = {};
end

switch command
   case 'on'
      if nargin ==1 || strcmp(fstr, 'all')
         eidors_cache('debug_on');
      end
      if ~eidors_debug('query',fstr);
         eidors_cache('debug_on', fstr);
      end
   case 'off'
      if nargin==1 || strcmp(fstr, 'all')
         eidors_cache('debug_off');
      else
         idx = strcmp(eidors_objects.debug_enabled_on, fstr);
         try
           eidors_objects.debug_enabled_on(idx) = [];
         end
      end
   case 'query'
      idx = strcmp(eidors_objects.debug_enabled_on, fstr);
      out = any(idx);
   otherwise
      error('EIDORS:WrongInput',['First input to eidors_debug must be ', ...
         '''on'',''off'' or ''query''.']);
end

function test
if eidors_debug('query','eidors_debug_test1')
   eidors_msg('TEST1',2);
end

if eidors_debug('query','eidors_debug_test2')
   eidors_msg('TEST2',2);
end


function do_unit_test
eidors_debug on eidors_debug_test1
eidors_debug on eidors_debug_test2
test
eidors_debug off eidors_debug_test1
test
eidors_debug off all
test
% expected output:
% TEST1
% TEST2
% TEST2
