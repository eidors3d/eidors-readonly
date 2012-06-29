function [eitimages_out,s] = EITCalcComplianceImage(eitimages)
%EITCALCCOMPLIANCEIMAGE   Computes a single dynamic compliance image
%surrogate from DeltaZ/DeltaP for the block of data.
%
% DeltaZ (%change from reference image) DeltaP cmH2O
%
% Signature:
% function [eitimages_out,s] = EITCalcComplianceImage(eitimages)
%
% Input:
%
% eitimages struct      scalar      .structure
%                                   .subject_ID
%                                   .subject_type
%                                   .case_date
%                                   .file_name
%                                   .eit_device
%                                   .image_series_length
%                                   .frame_rate
%                                   .difference_image_series
%                                   .absolute_iamge_series
%                                   .image_mask
%
% Output:
%
% eitimages_out struct      scalar      as eitimages_in with new members
%                                       .complianceimage
%
% s             boolean     scalar      errors present: true
%
% Copyright C. Gomez-Laberge, December 2010.
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
        help('EITCalcComplianceImage')
        display('Error EITCalcComplianceImage: invalid arguments');
        error('Error EITCalcComplianceImage: aborting execution');
end

if ~isfield(eitimages,'deltap')
    % Error
    s = true;
    help('EITCalcComplianceImage')
    display('Error EITCalcComplianceImage: Delta pressure unknown');
    error('Error EITCalcComplianceImage: aborting execution');
end
if ~isfield(eitimages,'fEIT')
    eitimages = EITCalcFunctionalImage(eitimages);
end

mask = ~eitimages.lungROI;
lungimage = eitimages.meantidalimage;
lungimage(mask)=-Inf;
C = lungimage/eitimages.deltap;
infmask = isinf(C);
negmask = C<0;
C(negmask & ~infmask)=0;

% Update output and return
eitimages_out = eitimages;
eitimages_out.complianceimage = C;

% End of function
end %function