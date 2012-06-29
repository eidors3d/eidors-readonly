function [vdout,s] = FormatViasysData(vdin)
%FORMATVIASYSDATA   Removes zeros between stimulation pairs in the Viasys
%voltage data.
%   
% Signature:
% function [vdout,s] = FormatViasysData(vdin)
%
% Input:
% vdin      double      MxN         Viasys voltage data with zeros
% 
% Output:   
% vdout     double      MxN         Viasys voltage data without zeros
% s         boolean     scalar      errors present: true
% 
% Copyright C. Gomez-Laberge, November 2010.
% $Id$

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
        help('FormatViasysData')
        display('Error FormatViasysData: invalid arguments');
        error('Error FormatViasysData: aborting execution');
end

indices = find(vdin(:,1)==0);
vdout = vdin;
checksum = sum(vdout(indices,:));
if checksum ~= 0
    % Report error without abort
    display('Error FunctionName: indices for zeros vary across frames--check output.');
    s = true;
end
vdout(indices,:) = [];

% End of function
end %function