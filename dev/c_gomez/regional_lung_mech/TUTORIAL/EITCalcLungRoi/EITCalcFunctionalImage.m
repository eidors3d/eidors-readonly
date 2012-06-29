function [eitimages_out,s] = EITCalcFunctionalImage(eitimages,range)
%EITCALCFUNCTIONALIMAGE   Computes a single functional image
%to the original image data. Functional images will only be computed from
%images within the specified range.
%
%The functional image is the pixel intensity standard deviation.
%
% Signature:
% function [eitimages_out,s] = EITCalcFunctionalImage(eitimages,range)
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
% (range)   double          1x2     [start,end] one segment only
%
% Output:
%
% eitimages_out struct      scalar      as eitimages_in with new members
%                                       .fEIT
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
        % Set (range) to default value
        range = [1,size(eitimages.difference_image_series,3)];
        option = 'none';
    case 2
        % Do nothing
    otherwise
        % Error
        s = true;
        help('EITCalcFunctionalImage')
        display('Error EITCalcFunctionalImage: invalid arguments');
        error('Error EITCalcFunctionalImage: aborting execution');
end

images = eitimages.difference_image_series;
images = images(:,:,range(1):range(2));
fEIT = std(images,0,3);

% Update output and return
eitimages_out = eitimages;
eitimages_out.fEIT = fEIT;

% End of function
end %function