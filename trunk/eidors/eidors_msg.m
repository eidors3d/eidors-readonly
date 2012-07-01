function log_level= eidors_msg( message, varargin )
%EIDORS_MSG: eidors progress and status messages
%
% USAGE: eidors_msg('mission accomplished; time to relax', 1)
%   prints 'EIDORS:[ mission accomplished; time to relax ]'
%       if log_level <=1
%   less important messages have higher levels
%
% other parameters will be sent to print; example
% eidors_msg('did %d of %s at %f', 2, 'stuff', sqrt(2), 1)
%   prints 'EIDORS:[ did 2 of stuff at 1.414214 ]'
%
% if the message starts with '@@@' then the character
%   is replaced with the calling function name
%
% USAGE: eidors_msg( 'log_level', 1)
%   sets the loglevel to 1
% USAGE: eidors_msg( 'log_level')
%   return the current log_level
%
% meanings of levels
%   0 => keep quiet (no messages should be level=0)
%   1 => important messages
%   2 => most messages
%   3 => detailed information
%
% Messages at log_level are displayed

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(message) && strcmp(message,'UNIT_TEST'); do_unit_test; return; end

global eidors_objects

if nargin==1
   args= {};
   level= 2;
else
   level= varargin{ nargin-1 };
   args= varargin( 1:nargin-2 );
end

[log_level] = get_levels;

for i= 1:length(args)
   if isa( args{i}, 'function_handle')
      args{i} = func2str(args{i});
   end
end

% It makes sense to print to stderr (fid=2), but matlab>7 prints this in red
ver= eidors_obj('interpreter_version');
if ver.ver>=7; fid= 1;
else ;         fid= 2; end


%try, lms= eidors_objects.last_message_size;
%  if lms>0; fprintf(fid,'%c', 8*ones(lms,1)); end
%end

% deal with variables
if ~ischar(message)
   var = inputname(1);
   message = [var ' = ' num2str(message)];
end


% Need to do twice to interpret text in message
% message= sprintf(message, args{:} );
if length(message)>1 % single characters are just for progress
   if length(message)>=3 && strcmp(message(1:3),'@@@');
      if length(message)==3; msg_extra = '';
      else                   msg_extra = [':', message(4:end)];
      end
      dbs = dbstack;
      message = sprintf('%s/%s%s', dbs(2).file, ...
            dbs(2).name, msg_extra);
   end
   string= [sprintf('%c',' ' * ones(1,level-1)), ...
            'EIDORS:[',message,']\n'];
       
elseif strcmp( message, 'log_level');
   log_level= eidors_objects.log_level;
else
   string= message;
end

if strcmp(message,'log_level')
   eidors_objects.log_level= level;
elseif level <= log_level
   fprintf(fid, string, args{:});
   if exist('OCTAVE_VERSION'); fflush(fid); end
   eidors_objects.last_message_size= 0;
end


function [log_level] = get_levels;
   global eidors_objects
   try
      log_level= eidors_objects.log_level;
   catch
      log_level= 2; % default;
      eidors_objects.log_level= log_level;
   end

function do_unit_test
   ll= eidors_msg('log_level');
   eidors_msg('log_level',5);
   eidors_msg('l1',1); eidors_msg('l2',2); eidors_msg('l3',3); eidors_msg('l4',4);

   eidors_msg('log_level',1);
   eidors_msg('l1',1); eidors_msg('l2',2); eidors_msg('l3',3); eidors_msg('l4',4);


   eidors_msg('@@@',1);
   eidors_msg('@@@ a message',1);
   eidors_msg('a @@@ message',1);
   extra_caller;

   eidors_msg('log_level',ll);

function extra_caller
   eidors_msg('@@@ a message from extra_caller',1);
  
