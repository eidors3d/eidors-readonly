function [vref,s] = CalcReferenceImage(vd,method)
%FUNCTIONTEMPLATE   Describe function here.
%   
% Signature:
% function [out1,out2] = FunctionTemplate(in1,in2,in3)
%
% Input:
% in1       double      scalar      description here
% in2       double      MxN         description here
% (in3)     boolean     scalar      description here
% 
% Output:   
% out1      double      scalar      description here
% out2      boolean     scalar      description here
% s         boolean     scalar      errors present: true
% 
% Copyright C. Gomez-Laberge, Month yyyy.
% $Id: $

% Set error status to 'no errors present'
s = false;
% Save calling path
callpath = pwd;

% Check arguments
switch nargin
    case 2
        % Do nothing
    otherwise
        % Error
        s = true;
        help('CalcReferenceImage')
        display('Error CalcReferenceImage: invalid arguments');
        error('Error CalcReferenceImage: aborting execution');
end

% % Define function constants
% K1=0;
% K2=1+i;

switch method
    case 'avgim'
        vref = mean(vd,2);
    case 'minim'
        vdsum = sum(vd);
        vref = vd(:,vdsum==min(vdsum));
    case 'maxim'
        vdsum = sum(vd);
        vref = vd(:,vdsum==max(vdsum));
    case 'leqavgim'
        vsum = sum(vd);
        vavgsum = mean(vsum);
        vd(:,vsum>vavgsum) = [];
        vref = mean(vd,2);
    case 'geqavgim'
        vsum = sum(vd);
        vavgsum = mean(vsum);
        vd(:,vsum<vavgsum) = [];
        vref = mean(vd,2);
    case 'minpixel'
        vref = min(vd')';
    case 'maxpixel'
        vref = max(vd')';
    case 'first10'
        vref = mean(vd(:,1:10),2);
    otherwise
        % Error
        s = true;
        help('CalcReferenceImage')
        display('Error CalcReferenceImage: invalid calculation method');
        error('Error CalcReferenceImage: aborting execution');
end

% % Report error with abort
% display('Error FunctionTemplate: error message');
% s = true;
% error('Error FunctionTemplate: aborting execution');
% 
% % Report error without abort
% display('Error FunctionName: error message');
% s = true;

% End of function
end %function