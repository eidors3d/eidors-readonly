function [eitimages_out,s] = EITCalcLungState(eitimages,Cmax,Pstar)
%EITCALCLUNGSTATE   Calculates proportion of regional alveolar collapse and
%over-distension relative to Cmax and Pstar, as described in 
%Costa et al., Intensive Care Med. (2009), 35:1132-1137.
%   
% Signature:
% function [eitimages_out,s] = EITCalcLungState(eitimages,Cmax,Pstar)
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
% Copyright C. Gomez-Laberge, December 2010.
% $Id$

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
        help('EITCalcLungState')
        display('Error EITCalcLungState: invalid arguments');
        error('Error EITCalcLungState: aborting execution');
end

% Define function constants
K1 = 10; % K1 results in percentile estimates with K1s precision

C = eitimages.complianceimage;

% Tabulate data using both NaN and lungROI
[Ctab,index] = TabulateImageData(C,eitimages.lungROI);
Cmaxtab = TabulateImageData(Cmax,eitimages.lungROI);
Pstartab = TabulateImageData(Pstar,eitimages.lungROI);
collapsetab = zeros(size(Ctab));
overdisttab = collapsetab;

% Calculate pixel maps
npixels = size(Ctab,1);
p = eitimages.peep;
for i = 1:npixels
    Cmaxtabi = Cmaxtab(i);
    Ctabi = Ctab(i);
    if Cmaxtabi <= 0 || Ctabi < 0
        collapsetab(i) = 0;
        overdisttab(i) = 0;
        continue
    end        
    if p < Pstartab(i)
        collapsetab(i) = (Cmaxtabi-Ctabi)/Cmaxtabi;
        overdisttab(i) = 0;
    elseif p > Pstartab(i)
        collapsetab(i) = 0;
        overdisttab(i) = (Cmaxtabi-Ctabi)/Cmaxtabi;
    else
        collapsetab(i) = 0;
        overdisttab(i) = 0;
    end
end
collapse = ImageTableData(collapsetab,size(C),index);
collapse(~eitimages.lungROI)=-Inf;
overdist = ImageTableData(overdisttab,size(C),index);
overdist(~eitimages.lungROI)=-Inf;

eitimages_out = eitimages;
eitimages_out.collapse = collapse;
eitimages_out.overdist = overdist;

% Calculate summary values
collapseSum = round(sum(collapsetab.*Cmaxtab)/sum(Cmaxtab)*100*K1)/K1;
overdistSum = round(sum(overdisttab.*Cmaxtab)/sum(Cmaxtab)*100*K1)/K1;
eitimages_out.collapseSum = collapseSum;
eitimages_out.overdistSum = overdistSum;

% End of function
end %function