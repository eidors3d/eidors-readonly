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
%
% meanings of levels
%   0 => keep quiet (no messages should be level=0)
%   1 => important messages
%   2 => most messages
%   3 => detailed information

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: eidors_msg.m,v 1.19 2007-08-29 09:25:00 aadler Exp $

global eidors_objects

if nargin==1
   args= {};
   level= 2;
else
   level= varargin{ nargin-1 };
   args= varargin( 1:nargin-2 );
end

try
   log_level= eidors_objects.log_level;
catch
   log_level= 2; % default;
   eidors_objects.log_level= log_level;
end

fid= 2; %stderr
if strcmp(message,'log_level')
   eidors_objects.log_level= level;
elseif level <= log_level
   string= sprintf('EIDORS:[ %s ]\n', message );
   fprintf(fid, string, args{:});
   if exist('OCTAVE_VERSION');
      fflush(fid);
   end
end
