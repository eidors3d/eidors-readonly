function [fpeak,fharmonic,s] = EITCalcFrequencySpectrum(eitdata,range,graph)
%EITCALCFREQUENCYSPECTRUM   Calculates and graphs the frequency spectrum of
%the voltage data averaged over all frames. Returns the frequency where the
%spectrum peak occurs.
%   
% Signature:
% function [fpeak,fharmonic,s] = EITCalcFrequencySpectrum(eitdata,range,graph)
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
% (range)   double      1x2         [start,end] one segment only
% (graph)   boolean     scalar      graph FFT if true
% 
% Output:   
% fpeak     double      scalar      frequency (Hz) of spectral peak
% fharmonic double      scalar      distance to next harmonic
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
        % Set (range) to default value
        range = [1,size(eitdata.voltage_data,2)];
        % Set (graph) to default value
        graph = true;
    case 2
        % Set (graph) to default value
        graph = true;
    case 3
        % Do nothing
    otherwise
        % Error
        s = true;
        help('EITCalcFrequencySpectrum')
        display('Error EITCalcFrequencySpectrum: invalid arguments');
        error('Error EITCalcFrequencySpectrum: aborting execution');
end

% % Define function constants
K1=0.15; % Lower cutoff in Hz

% Extract zero-mean voltage data
v0 = sum(eitdata.voltage_data(:,range(1):range(2)),1);
if isempty(eitdata.frame_rate)
        s = true;
        display('Error EITCalcFrequencySpectrum: frame rate unspecified');
        error('Error EITCalcFrequencySpectrum: aborting execution');
end    
[Ymag,Ypha,f] = CalcFourierTransform(v0,eitdata.frame_rate,graph);
lowercutoff = find(f>K1,1);
trunkYmag = Ymag;
trunkYmag(1:lowercutoff)=0;
fmax = find(trunkYmag == max(trunkYmag));
fpeak = f(fmax);
fharmonic = f(fmax-1);

% End of function
end %function