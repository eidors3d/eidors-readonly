function [vi,vq,s] = GetQuadratureSignals(vd)
%GETQUADRATURESIGNALS   Function extracts the two quadrature signals
%recorded in the Draeger EIT device.
%   
% Signature:
% function [vi,vq] = GetQuadratureSignals(vd)
%
% Input:
% vd        double      MxN      raw voltage data
% 
% Output:   
% vi        double      MxN      positive quadrature voltage signal
% vq        double      MxN      negative quadrature voltage signal
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
        help('GetQuadratureSignals')
        display('Error GetQuadratureSignals: invalid arguments');
        error('Error GetQuadratureSignals: aborting execution');
end

% Determine file version
KV32        = 513;   % frame length per version
KV401       = 645;
KV32ISTART  = 2;     % start and stop indices for vi and vq per version
KV32IEND    = 209;
KV32QSTART  = 258;
KV32QEND    = 465;
KV401ISTART  = 2;
KV401IEND    = 209;
KV401QSTART  = 324;
KV401QEND    = 531;

% Get signals according to file version
switch size(vd,1)
    case KV32
        vi = vd(KV32ISTART:KV32IEND,:);
        vq = vd(KV32QSTART:KV32QEND,:);
    case KV401
        vi = vd(KV401ISTART:KV401IEND,:);
        vq = vd(KV401QSTART:KV401QEND,:);
    otherwise
        s = true;
        display('Error GetQuadratureSignals: unknown file version');
        error('Error GetQuadratureSignals: aborting execution');
end

% End of function
end %function