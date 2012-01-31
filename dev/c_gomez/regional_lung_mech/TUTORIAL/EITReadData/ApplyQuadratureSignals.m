function [vd,s] = ApplyQuadratureSignals(vi,vq,code)
%APPLYQUADRATURESIGNALS   Computes output signal given quadrature signals
%according to code supplied:
%   'vi'  - the positive quadrature signal
%   'vq'  - the negative quadrature signal
%   'sum' - the sum of the quadrature signals
%   
% Signature:
% function vd = ApplyQuadratureSignals(vi,vq,code)
%
% Input:
% vi       double      MxN         the positive quadrature signal
% vq       double      MxN         the negative quadrature signal
% code     char        Nx1         the code word
% 
% Output:   
% vd        double      MxN         description here
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
    case 3
        % Do nothing
    otherwise
        % Error
        s = true;
        help('ApplyQuadratureSignals')
        display('Error ApplyQuadratureSignals: invalid arguments');
        error('Error ApplyQuadratureSignals: aborting execution');
end

% Compute output signal according to code supplied
switch code
    case 'vi'
        vd = vi;
    case 'vq'
        vd = vq;
    case 'sum'
        vd = vi+vq;
    otherwise
        % Error
        s = true;
        help('ApplyQuadratureSignals')
        display('Error ApplyQuadratureSignals: unknown code');
        error('Error ApplyQuadratureSignals: aborting execution');
end

% End of function
end %function