function [eitimages_out,optimalp,s] = EITCalcLungRoi(eitimages,optimalp)
%EITCALCLUNGROI   Determines lung ROI using fEIT images. Determines cutoff
%value for ROI heuristically: threshold perturbation causing least effect.
%
% Copyright C. Gomez-Laberge, December 2010.
% $Id: $

s = false;

switch nargin
    case 1
        optimalp=[];
    case 2
        % Do nothing
    otherwise
        % Error
        help('EXELungRoi')
        display('Error EXELungRoi: invalid arguments');
        error('Error EXELungRoi: aborting execution');
end

if ~isfield(eitimages,'fEIT')
    eitimages = EITCalcFunctionalImage(eitimages);
end

DEBUG = 0;

parray = [0:0.05:1];
feit = eitimages.fEIT;
sdmax = max(feit(:));

if isempty(optimalp)
    optimalerror = 1;
    for i = 2:length(parray)-1
        roineg = feit > parray(i-1)*sdmax;
        roitest = feit > parray(i)*sdmax;
        roipos = feit > parray(i+1)*sdmax;
        errorneg = (nnz(roineg)-nnz(roitest))/nnz(roitest);
        errorpos = (nnz(roitest)-nnz(roipos))/nnz(roipos);
        maxerror = max(errorneg,errorpos);
        if maxerror < optimalerror
            optimalerror = maxerror;
            optimalp = parray(i);
        end
        if DEBUG
            figure, imagesc(roitest), colormap(gray);
            title(['p=' num2str(parray(i))]);
            axis square off
        end
    end
end

lungROI = feit > optimalp*sdmax;
eitimages_out = eitimages;
eitimages_out.lungROI = lungROI;
% End of function
end %function