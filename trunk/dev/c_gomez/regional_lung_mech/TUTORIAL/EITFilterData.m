function [eitdata_out,s] = EITFilterData(eitdata_in,filtertype,range,option)
%EITFILTERDATA   Builds a filter of given type, filters voltage data, and
%returns a new eitdata structure with filtered voltage data and a filter
%member listing filter parameters used.
%
%
% lowpass to remove cardiac signal (a zero-phase FIR equiripple filter)
% bandpass to remove baseline drift and cardiac (not yet implemented)
% highpass to remove baseline and breathing (not yet implemented
%
% Signature:
% function [eitdata_out,s] = EITFilterData(eitdata_in,filtertype,range,option)
%
% Input:
% eitdata_in    double      struct      unfiltered eitdata
% filtertype    char        1xN         {'lowpass','bandpass','highpass',
%                                        'ad hoc'}
% (range)       double      1x2         [start,end] one segment only
% (option)      char        1xN         {'none','manual'}
%
% Output:
% eitdata_out   double      struct      filtered eitdata with new member
%                                       .filterHd (visualize with fvtool)
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
        range = [1,size(eitdata_in.voltage_data,2)];
        option = 'none';
    case 3
        % Set (range) to default value
        option = 'none';
    case 4
        % Do nothing
    otherwise
        % Error
        s = true;
        help('EITFilterData')
        display('Error EITFilterData: invalid arguments');
        error('Error EITFilterData: aborting execution');
end

% Define function constants
K1=100;  % filter order
K2=1;   % number of harmonics after peak for transition band
K3 = 0.6641; % Manual frequency peak setting
K4 = 0.9961-K3; % Manual frequency harmonic setting
%     K3 = 0.6094;    % 1006-RERV c1
%     K4 = 1.1930 -K3; % 1006-RERV c1
%     K3 = 0.4063;    % 1004-RERV d5
%     K4 = .8379 -K3; % 1004-RERV d5
K5 = 0.1; % Bandpass Fstart/fpeak ratio

% The sampling rate
Fs = eitdata_in.frame_rate;

% Determine filter type
switch filtertype
    case 'lowpass'
        % Find spectral peak
        if strcmp(option,'manual')
            fpeak = K3;
            fharmonic = K4;
        else
            [fpeak,fharmonic,s] = EITCalcFrequencySpectrum(eitdata_in,range,false);
        end
        % Calculate transition band as K2+1 harmonic of peak
        knee = fpeak + K2*fharmonic;
        d = fdesign.lowpass('N,Fc',K1,knee,Fs);
        filterHd = design(d);
        % Filter data
        fvd = filtfilt(filterHd.Numerator,1,eitdata_in.voltage_data')';
    case 'bandpass'
        % Find spectral peak
        if strcmp(option,'manual')
            fpeak = K3;
            fharmonic = K4;
        else
            [fpeak,fharmonic,s] = EITCalcFrequencySpectrum(eitdata_in,range,false);
        end
        Fstop1 = 0.01; %fpeak*K5;
        Fstop2 = fpeak + K2*fharmonic;
        % Input of parameters
        d = fdesign.bandpass('n,fc1,fc2');
        d = fdesign.bandpass('n,fc1,fc2',K1,Fstop1,Fstop2,Fs);        
        filterHd = design(d);
        % Filter data
        fvd = filtfilt(filterHd.Numerator,1,eitdata_in.voltage_data')';
    case 'highpass'
        display('Error EITFilterData: highpass not implemented');
        s = true;
        error('Error EITFilterData: aborting execution');
    case 'ad hoc'
        fvd = find_artefacts(eitdata_in.voltage_data);
    otherwise
        display('Error EITFilterData: unknown filter type');
        s = true;
        error('Error EITFilterData: aborting execution');
end

% Create output structure
eitdata_out = eitdata_in;
eitdata_out.voltage_data = fvd;
if ~strcmp(filtertype,'ad hoc')
    eitdata_out.filterHd = filterHd;
end

% End of function
end %function