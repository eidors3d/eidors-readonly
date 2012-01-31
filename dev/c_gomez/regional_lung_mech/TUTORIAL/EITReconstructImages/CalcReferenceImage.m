function [vref,s] = CalcReferenceImage(vd,method)
% Copyright C. Gomez-Laberge, January 2012.
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

% End of function
end %function