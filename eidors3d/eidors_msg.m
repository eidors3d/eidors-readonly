function eidors_msg( message, varargin )
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
% USAGE: eidors_msg( 'log_level', 1)
%   sets the loglevel to 1
% USAGE: eidors_msg( 'log_level2', 3)
%   sets the loglevel2 to 3
%
% meanings of levels
%   0 => keep quiet (no messages should be level=0)
%   1 => important messages
%   2 => most messages
%   3 => detailed information
%
% Messages at log_level are displayed
% Messages at log_level2 are displayed temprorarily

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: eidors_msg.m,v 1.26 2008-03-28 19:55:38 aadler Exp $

global eidors_objects

if nargin==1
   args= {};
   level= 2;
else
   level= varargin{ nargin-1 };
   args= varargin( 1:nargin-2 );
end

[log_level, log_level2] = get_levels;

for i= 1:length(args)
   if isa( args{i}, 'function_handle')
      args{i} = func2str(args{i});
   end
end

fid= 2; %stderr
%try, lms= eidors_objects.last_message_size;
%  if lms>0; fprintf(fid,'%c', 8*ones(lms,1)); end
%end

% Need to do twice to interpret text in message
message= sprintf(message, args{:} );
if length(message)>1 % single characters are just for progress
   string= sprintf('EIDORS:[ %s ]\n', message);
else
   string= message;
end

if strcmp(message,'log_level')
   eidors_objects.log_level= level;
elseif level <= log_level
   fprintf(fid, string);
   if exist('OCTAVE_VERSION'); fflush(fid); end
   eidors_objects.last_message_size= 0;
elseif 0 %level <= log_level2
% TODO, save messages, and write them out via (printing-hyperlinks)
%   so that users can see them.
   fprintf(fid, string);
   eidors_objects.last_message_size= length(string);
end


function [log_level, log_level2] = get_levels;
   global eidors_objects
   try
      log_level= eidors_objects.log_level;
   catch
      log_level= 2; % default;
      eidors_objects.log_level= log_level;
   end
   try
      log_level2= eidors_objects.log_level2;
   catch
      log_level2= 4; % default;
      eidors_objects.log_level2= log_level2;
   end

