function [status,result] = system_cmd( cmd );
% SYSTEM_CMD: issue system commands, and try to compensate for
%    strange differences between systems.
%
% status = system_cmd( cmd );

% (C) 2012 Andy Adler, License GPL v2 or v3
% $Id$

if ~exist('OCTAVE_VERSION')
   if strfind(system_dependent('getos'),'Linux')
     %Version 7 under linux sets the LD_LIBRARY_PATH and that breaks netgen    
     ldpath ='LD_LIBRARY_PATH=;';
     % Problems are a) some newer matlab versions dont' need it,
     %   b) what if people are working in tcsh?

     cmd = ['sh -c ''LD_LIBRARY_PATH=""; ',cmd,' '' '];
   end
end

[status,result] = system(cmd);
