function eidors_msg( message, level )
%EIDORS_MSG: eidors progress and status messages
%
% USAGE: eidors_msg('mission accomplished; time to relax', 1)
%   prints 'mission accomplished; time to relax' if loglevel <=1
%   less important messages have higher levels
%
% USAGE: eidors_msg( 'set_level', 1)
%   sets the loglevel to 1
%
% meanings of levels
%   0 => keep quiet (no messages should be level=0)
%   1 => important messages
%   2 => most messages
%   3 => detailed information

global eidors_objects

try
   log_level= eidors_objects.log_level;
catch
   log_level= 1; % default;
   eidors_msg.log_level= log_level;
end

fid= 2; %stderr
if level <= log_level
   fprintf(fid,  message );
   if exist('OCTAVE_VERSION');
      fflush(fid);
   end
end
