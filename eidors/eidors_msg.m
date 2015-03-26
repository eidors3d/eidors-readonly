function log_level= eidors_msg( message, varargin )
%EIDORS_MSG eidors progress and status messages
%
%  EIDORS_MSG(STR, LVL) prints 'EIDORS:[ STR ]' if current log_level <= LVL
%     If the message starts with '@@@' then the characters are replaced 
%     with the calling function's name. Default LVL = 2.
%
%  EIDORS_MSG(STR, LVL, ...) uses STR as format string for fprintf passing
%     to it any additional arguments. 
%
%  LVL = EIDORS_MSG('log_level') returns the current log_level
%
%  EIDORS_MSG('log_level', LVL) sets the log_level for all of EIDORS
%
%  EIDORS_MSG('log_level', LVL, FNAME) sets a custom log_level for m-file
%     FNAME (provide without extension)
%
%  EIDORS_MSG('log_level','reset', FNAME) resets the log_level for m-file
%     FNAME to the default level 
%
%  EIDORS_MSG('log_level','reset')
%  EIDORS_MSG('log_level','reset','all') resets all custom log-level
%     settings to default
%
%  Meanings of levels:
%     0 => keep quiet (no messages should be level=0)
%     1 => important messages
%     2 => most messages
%     3 => detailed information%
%
%  Examples:
%  eidors_msg('mission accomplished; time to relax', 1) prints 
%     >> EIDORS:[ mission accomplished; time to relax ]'
%       if log_level <=1
%  less important messages have higher levels
%
%  eidors_msg('did %d of %s at %f', 2, 'stuff', sqrt(2), 1)
%     >> EIDORS:[ did 2 of stuff at 1.414214 ]
% 
%  eidors_msg('@@');
%     >> EIDORS:[eidors_msg]
%
%  eidors_msg('@@@');
%     >> EIDORS:[eidors_msg/do_unit_test]
%
%  eidors_msg('@@@ a message',1);
%     >> EIDORS:[eidors_msg/do_unit_test: a message]
%
%  See also EIDORS_DEBUG, EIDORS_CACHE.

% (C) 2005-2014 Andy Adler and Bartlomiej Grychtol
% License: GPL version 2 or version 3
% $Id$

if ischar(message) && strcmp(message,'UNIT_TEST'); do_unit_test; return; end

if ischar(message) && strcmp(message, 'log_level')
   log_level = process_log_level(varargin{:});
   return;
end

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
   if length(message)>=2 && strcmp(message(1:2),'@@');
      if length(message)==2; msg_extra = '';
      else                   msg_extra = [':', message(3:end)];
      end
      [file fun] = get_caller();
      dbs = dbstack;
      if length(message) > 2 && message(3) == '@'
          if length(message) == 3
              msg_extra = '';
          else
              msg_extra(2) = [];
          end
          message = sprintf('%s/%s%s', file, fun, msg_extra);
      else
          message = sprintf('%s%s', file , msg_extra);
      end
   end
   string= [sprintf('%c',' ' * ones(1,level-1)), ...
            'EIDORS:[',message,']\n'];
       
else
   string= message;
end

if level <= log_level
   fprintf(fid, string, args{:});
   if exist('OCTAVE_VERSION'); fflush(fid); end
   eidors_objects.last_message_size= 0;
end


function log_level = process_log_level(varargin)
   global eidors_objects
   if isempty(varargin)
      % eidors_msg('log_level');
      log_level = eidors_objects.log_level;
      return
   end
   if isnumeric(varargin{1})
      log_level = eidors_objects.log_level;
      switch nargin
         case 1
            % eidors_msg('log_level', 2);
            eidors_objects.log_level = varargin{1};
         case 2
            % eidors_msg('log_level', 1, fname);
            caller = get_caller;
            custom = eidors_objects.log_level_custom;
            if isempty(custom), return, end
            idx = find(ismember(caller, custom.names), 1);
            if isempty(idx), idx = length(custom.levels) + 1; end
            eidors_objects.log_level_custom.names{idx} = caller;
            eidors_objects.log_level_custom.levels(idx) = varargin{1};
         otherwise
            error('Wrong number of inputs');
         
      end
   elseif ischar(varargin{1}) && strcmp(varargin{1}, 'reset')
      log_level = eidors_objects.log_level;
         
      switch nargin
         case 1
            % eidors_msg('log_level','reset');
            eidors_objects.log_level_custom.names = {};
            eidors_objects.log_level_custom.levels = []; 
         case 2
            switch varargin{2}
               case 'all'
                  % eidors_msg('log_level','reset', 'all');
                  eidors_objects.log_level_custom.names = {};
                  eidors_objects.log_level_custom.levels = [];
               otherwise
                  caller = get_caller;
                  idx = find(ismember(caller, ...
                              eidors_objects.log_level_custom.names{1}), 1);
                  eidors_objects.log_level_custom.names(idx) = [];
                  eidors_objects.log_level_custom.levels(idx) = [];
            end
         otherwise
            error('Wrong number of inputs');
      end
            
   else
      error('Wrong input, consult help.');
      
   end
      

function [file fun] = get_caller
   dbs = dbstack;
   file = dbs(3).file;
   idx = find(file == '.', 1, 'last');
   if isempty(idx), idx = length(file)+1; end
   file = file(1:idx-1);

   fun =  dbs(3).name;

function [log_level] = get_levels
   global eidors_objects
   try 
      custom = eidors_objects.log_level_custom;
   catch
      custom.names = {};
      custom.levels = [];
      eidors_objects.log_level_custom = custom;
   end
   if ~isempty(custom)
      % find caller file name
      file = get_caller;
      % check for custom flag
      idx = find(ismember(file, custom.names));
      if ~isempty(idx)
         log_level = custom.levels(idx);
         return
      end
   end
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

   eidors_msg('@@',1);
   eidors_msg('@@@',1);
   eidors_msg('@@ a message',1);
   eidors_msg('@@@ a message',1);
   eidors_msg('a @@ message',1);
   eidors_msg('a @@@ message',1);
   extra_caller;

   % test custom levels
   eidors_msg('log_level',2);
   eidors_msg('FAIL',3);
   eidors_msg('log_level',1, 'eidors_msg');
   eidors_msg('PASS',1);
   eidors_msg('FAIL',2);
   eidors_msg('log_level',3, 'eidors_msg');
   eidors_msg('PASS',3);
   eidors_msg('FAIL',4);
   eidors_msg('log_level','reset','eidors_msg');
   eidors_msg('PASS',2);
   eidors_msg('FAIL',3);
   % just for fun
   eidors_msg('log_level','reset'); 
   eidors_msg('log_level','reset','all'); 
   
   eidors_msg('log_level',ll);
   

function extra_caller
   eidors_msg('@@@ a message from extra_caller',1);
   eidors_msg('@@ a shorter message from extra caller',1);
  
