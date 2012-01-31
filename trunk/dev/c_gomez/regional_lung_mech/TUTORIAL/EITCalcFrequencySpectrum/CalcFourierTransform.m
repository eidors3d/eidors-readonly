function [Ymag,Ypha,f] = CalcFourierTransform(y,samplerate,graph)
%CALCFOURIERTRANSFORM   Calculates and graphs the FFT of signal y.
%
% Signature:
% function [Ymag,Ypha] = CalcFourierTransform(y,samplerate,graph)
%
% Input:
% y             double      1xN         signal
% samplerate    double      scalar      sampling rate
% (graph)       boolean     scalar      graph FFT if true
%
% Output:
% Ymag          double     1xM          magnitude of FFT
% Ypha          double     1xM          phase of FFT
% f             double     1xM          frequency axis (Hz)
%
% Copyright C. Gomez-Laberge, November 2010.
% $Id: $

% Check arguments
switch nargin
    case 2
        % Set (graph) to default value
        graph = true;
    case 3
        % Do nothing
    otherwise
        % Error
        help('CalcFourierTransform')
        display('Error CalcFourierTransform: invalid arguments');
        error('Error CalcFourierTransform: aborting execution');
end

% Define function constants
% K1=0;
% K2=1+i;

% Remove signal offset
y = y-mean(y);
L = length(y);
t = (0:L-1)/samplerate;
nfft = 2^nextpow2(L);
Y = fft(y,nfft)/L;
f = samplerate/2*linspace(0,1,nfft/2+1);
Ymag = 2*abs(Y(1:nfft/2+1));
Ypha = 2*angle(Y(1:nfft/2+1));

if graph
    figure, subplot(2,1,1), plot(f,Ymag)
    ylabel('|Y(f)|')
    subplot(2,1,2), plot(f,Ypha)
    ylabel('arg Y(f)')
    xlabel('frequency (Hz)')
end

% End of function
end %function