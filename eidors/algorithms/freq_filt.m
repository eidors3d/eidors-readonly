function vfilt = freq_filt(vin, fresp, FR, dim)
%FREQ_FILT: frequency filter data
% vfilt = freq_filt(vin, fresp, FR, dim)
% vin = matrix of data values
%     = data structure (with field vin.meas)
%     = image structure (with field vin.elem_data)
% fresp = function of freq filter
%      e.g. fresp = @(f) f<10;  % values in Hz
% FR  = frame rate (or use FR parameter on data structures)
% dim = dimension along which to filter (default is 2)

% (C) Andy Adler 2019. License: GPL v2 or v3.
% $Id$

% if FR not given, see if it's a parameter

% TODO: Add windowing

if ischar(vin) && strcmp(vin,'UNIT_TEST'); do_unit_test; return; end

if nargin<3; 
   FR = vin.FR;
end
if nargin<4; 
   dim = 2;
end

if isstruct(vin) && isfield(vin,'type');
   switch vin.type
     case 'data'
        vfilt.meas = do_freq_filt(vin.meas, fresp, FR, dim);
     case 'image'
        vfilt.elem_data = do_freq_filt(vin.elem_data, fresp, FR, dim);
     otherwise
        error('Can''t process object of type %',vin.type);
   end
elseif isnumeric(vin)
   vfilt = do_freq_filt(vin, fresp, FR, dim);
else 
   error('can''t process object');
end
      


% Filter in the frequency direction
function s = do_freq_filt(s,fresp,FR, dim);
  f = fft(s,[],dim);
  fshape = [1,1,1]; fshape(dim) = size(s,dim);
  fax = freq_axis_filt(FR,size(s,dim),fresp);
  f = f .* reshape(fax,fshape);
  s= ifft(f,[],dim);
  if norm(imag(s(:))) > 1e-11
     error('FFT filter has imag output');
  end
  s = real(s);
  
function fax = freq_axis_filt( FR, lD, fresp);
  fax = linspace(0,FR,lD+1);
  fax(end)=[];
  fax(fax>FR/2) = fax(fax>FR/2) - FR;
  fax = feval(fresp, abs(fax));

function do_unit_test
  FR = 100; t = (0:1e3)/FR;
  s = sin(2*pi*5*t) + 2*cos(2*pi*15*t) +  3*sin(2*pi*0.1*t);
  plot(t,s); hold on;
 
  for fnum = 1:3; switch fnum
     case 1; fresp = @(f) f<10;
     case 2; fresp = @(f) f<1;
     case 3; fresp = @(f) (f<1) + (f>=1)./(f+eps);
     end
     sf= freq_filt(s,fresp, FR);
     plot(t,sf,'LineWidth',2);
  end
  hold off;
