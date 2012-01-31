function [eitdata_out,eitimages_out,s] = EITCalcTidalImages(eitdata,eitimages,range,option)
%EITCALCTIDALIMAGES   Computes a series of tidal images and their indices
%to the original image data. If a time range is
%specified, then tidal images will only be constructed within. The output
%is a copy of the input with two new members:
%.tidalimages and .tidalindices.
%
%A tidal breath is computed by subtracting the most recent image
%corresponding to a local minima from the current image at a local maxima.
%
%Option = 'flip' swaps minima-maxima times
%
% Signature:
% function [eitdata_out,eitimages_out,s] = EITCalcTidalImages(eitdata,eitimages,range,option)
%
% Input:
% eitdata   struct      scalar      .structure
%                                   .subject_ID
%                                   .subject_type
%                                   .case_date
%                                   .file_name
%                                   .eit_device
%                                   .number_of_frames
%                                   .frame_rate
%                                   .voltage_data
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
% (option)  char            1xN     {'none','flip'}
%
% Output:
% eitdata_out   struct      scalar      .structure
%                                       .subject_ID
%                                       .subject_type
%                                       .case_date
%                                       .file_name
%                                       .eit_device
%                                       .number_of_frames
%                                       .frame_rate
%                                       .voltage_data
%                                       .tidalindices
%
% eitimages_out struct      scalar      as eitimages_in with new members
%                                       .tidalimages
%                                       .tidalindices
%
% s             boolean     scalar      errors present: true
%
% Copyright C. Gomez-Laberge, November 2010.
% $Id: $

% Set error status to 'no errors present'
s = false;
% Save calling path
callpath = pwd;

% Check arguments
switch nargin
    case 2
        % Set (range) to default value
        range = [1,size(eitdata.voltage_data,2)];
        option = 'none';
    case 3
        option = 'none';
    case 4
        % Do nothing
    otherwise
        % Error
        s = true;
        help('EITCalcTidalImages')
        display('Error EITCalcTidalImages: invalid arguments');
        error('Error EITCalcTidalImages: aborting execution');
end

% Define function constants
GLOBALSIGNAL = 1; % set to 1 for voltage data, 0 for image data
K1=0.5; % tidal period fraction for tidalperiodthreshold
K2 = 1; % tidal SD fraction for tidalmagnitudethreshold

% Calc global signal nu
if GLOBALSIGNAL
    nufull = mean(eitdata.voltage_data,1);
else
    nufull = EITCalcTimeSignal(eitimages,[],0);
end    
if strcmp(option,'flip')
    ntemp = -nufull;
    ntemp = ntemp + 2*mean(nufull(:));
    nufull = ntemp;
end

% Select segments in range
nu = nufull(range(1):range(2));
% Differentiate twice to find local extrema
Dnu = [0 diff(nu)];
signDnu = sign(Dnu)<0;
DsignDnu = [0, diff(signDnu)];
% Discard partial breaths at start and end if present
maxima = find(DsignDnu == 1);
minima = find(DsignDnu == -1);
fpeak = EITCalcFrequencySpectrum(eitdata,range,false);
tidalperiodthreshold = K1*eitdata.frame_rate/fpeak;
if maxima(1) < tidalperiodthreshold
    maxima(1) = [];
elseif minima(1) > maxima(1)
    minima = [1,minima];
end
if length(nu)-maxima(end) > tidalperiodthreshold && maxima(end) < minima(end)
    minima(end) = [];
end

% Readjust extema value to match original signal nufull
maxima = maxima+range(1)-1;
minima = minima+range(1)-1;

% Step through each non-zero value for the extrema voltages
% and decide if adjacent vmax - vmin is a true tidal breath based on the
% nu signal difference at those points EG (if vmax-vmin < std(nu)) then
% ignore the pair
tidalmagnitudethreshold = K2*std(nu);
for i = 1:length(maxima)
    if nufull(maxima(i))-nufull(minima(i)) < tidalmagnitudethreshold
        maxima(i) = NaN;
    end
end
deleteindices = isnan(maxima);
maxima(deleteindices) = [];
minima(deleteindices) = [];
if length(minima)>length(maxima)
    minima(end)=[];
end

% Calculate all tidal images
breaths = length(maxima);
zdelta = eitimages.difference_image_series;
tidalimagesize = size(zdelta);
tidalimagesize(3) = breaths;
tidalimages = zeros(tidalimagesize);
for i = 1:breaths
    temp = zdelta(:,:,maxima(i)) - zdelta(:,:,minima(i));
    tidalimages(:,:,i) = temp;
end

% Update output and return
eitdata_out = eitdata;
eitdata_out.tidalindices = [maxima;minima];
eitimages_out = eitimages;
eitimages_out.tidalindices = [maxima;minima];
eitimages_out.tidalimages = tidalimages;
eitimages_out.meantidalimage = mean(tidalimages,3);

% End of function
end %function