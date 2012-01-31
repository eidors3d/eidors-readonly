function [type,s] = GetSubjectType(studycode)
%GETSUBJECTTYPE   Cross references studycode with known studies and returns
%subject type 'Human', 'Swine', 'Simulation'.
%   
% Signature:
% function [type,s] = GetSubjectType(studycode)
%
% Input:
% studycode char        Nx1
% 
% Output:   
% type      char        Nx1         {'Human', 'Swine', 'Simulation', 
%                                    'Test', 'Unknown'}
% s         boolean     scalar      errors present: true
% 
% Copyright C. Gomez-Laberge, November 2010.
% $Id: $

% Set error status to 'no errors present'
s = false;
% Save calling path
callpath = pwd;

% Check arguments
switch nargin
    case 1
        % Do nothing
    otherwise
        % Error
        s = true;
        help('GetSubjectType')
        display('Error GetSubjectType: invalid arguments');
        error('Error GetSubjectType: aborting execution');
end

switch studycode
    case 'STUDYNAME'
        type = 'Human';
    case 'TEST'
        type = 'Test';   
    otherwise
        s = false;
        type = 'Unknown';
        display('Error GetSubjectType: unknown study code');
end
        
% End of function
end %function