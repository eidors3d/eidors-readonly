function [C,th] = fourier_fit(points,N,start);
% FOURIER_FIT: use fourier series to interpolate onto a boundary
%
% [pp] = fourier_fit(points) fits a Fourier series
%    points - [x y] contour to be fitted
% [pp] = fourier_fit(points,N) fits a Fourier series and downsamples
%    N is the number of Fourier components to fit to.
%
% [xy] = fourier_fit(pp,  linear_frac, start) returns points at linear_frac
% distance along the contour
%    pp          - piecewise polynomial structure
%    linear_frac - vector of length fractions (0 to 1) to calculate points
%    start       - (optional) a seed for the first point
%    xy          - cartesian coordinates of the points
%
%
% xy = fourier_fit(points,'fit_spacing',spacing);
% example
%  xy = fourier_fit(point,'fit_spacing',0:.1:999);
%  th = atan2(xy(:,2)-centroid(2),xy(:,1)-centroid(1));
% 

% (C) Andy Adler, 2010. Licenced under GPL v2 or v3
% $Id$

if ischar(points) && strcmp(points,'UNIT_TEST'); do_unit_test; return ; end

if nargin<2; N= size(points,1); end
if ischar(N) && strcmp(N,'fit_spacing')
   pp = fourier_fit(points, size(points,1),points(1,:));
   C = fourier_fit(pp, start);
   return
end


if size(points,2)==2 % calling analysis function
   C = fit_to_fourier(points,N);
elseif size(points,2)==1 % calling synthesis function
   if nargin<3; start = []; end
   C = fit_from_fourier(points,N,start);
else
   error('size of first parameter to fourier_fit not undersood');
end

% this will crash if N>length(points)
function C = fit_to_fourier(points,N);
   z = points*[1;1i];
   Z = fft(z,max(N,size(points,1)));
   if rem(N,2)==0 % Even 
     N2 = N/2;
     Zlp = Z([1,2:N2,1,end-(N2-2:-1:0)]);
     Zlp(N2+1) = 0.5*(Z(N2+1) + Z(end-N2+1)); %centre point
   else 
     N2 = floor(N/2);
     Zlp = Z([1,2:N2+1,end-(N2-1:-1:0)]);
   end
   C = length(Zlp)/length(Z)*Zlp; % Amplitude scaling
    
function xy = fit_from_fourier(C,linear_frac,start);
   % Step 1: oversample
   N = length(C);
   n2 = ceil(N/2);

   pad = zeros(10000,1);
   if rem(N,2)==0 % even
      Zos = [C(1:n2); C(n2+1)/2; pad; C(n2+1)/2; C(n2+2:end)];
   else
      Zos = [C(1:n2); pad; C(n2+1:end)];
   end
   Zos = length(Zos)/length(C)*Zos;
   zos = ifft(Zos);
   % rearrange the points such that they start as close as possible to the
   % seed
   if ~isempty(start)
       if size(start,2) == 1; start = start'; end
       start = start*[1; 1i];
       dist = abs(zos-start);
       [jnk,p] = min(dist);
       zos = circshift(zos,-p+1);
   end
   
   % Step 2:
   zos(end+1) = zos(1); % make sure the loop is closed
   dpath= abs(diff(zos));
   pathlen = [0;cumsum(dpath)];

   % interpolate points onto path
   npath = pathlen/max(pathlen);
   linear_frac= mod(linear_frac,1);
   z_int = interp1(npath, zos, linear_frac);

   xy= [real(z_int(:)), imag(z_int(:))];
   
function do_unit_test

a=[-0.8981  0.1404;-0.7492  0.5146;-0.2146  0.3504;
    0.3162  0.5069; 0.7935  0.2702; 0.9615 -0.2339;
    0.6751 -0.8677; 0.0565 -0.6997;-0.3635 -0.8563;
   -0.9745 -0.4668];

   C= fourier_fit(a,10);
   xy = fourier_fit(C, linspace(.05,1.04,100));
   xy2= fourier_fit(C, [.3,.4,.5,.6]);
   plot(a(:,1),a(:,2),'r',xy(:,1),xy(:,2),'b',xy2(:,1),xy2(:,2),'m*');

   a(5,:)= [];
   C= fourier_fit(a);
   xy = fourier_fit(C, linspace(.05,1.04,100));
   xy2= fourier_fit(C, [.3,.4,.5,.6]);
   plot(a(:,1),a(:,2),'r',xy(:,1),xy(:,2),'b',xy2(:,1),xy2(:,2),'m*');
   eidors_msg('VIEW GRAPH TO VERIFY',0);

subplot(212);
 n_elecs = 8; centroid = mean(a);
 p = linspace(0,1,n_elecs+1)'; p(end) = [];
 xy = fourier_fit(a,'fit_spacing',p);
 plot(a(:,1),a(:,2),'r',xy(:,1),xy(:,2),'b*');

  centroid = mean(a);
 th = atan2(xy(:,2)-centroid(2), xy(:,1)-centroid(1));
 hold on;
 aa = [xy(:,1),centroid(1)+0*xy(:,1)]';
 bb = [xy(:,2),centroid(1)+0*xy(:,2)]';
 plot(aa,bb,'k-');
  unit_test_cmp('fit th',diff(th(1:4)), [ -0.671324440246111;
  -1.011260720375461; -0.791975282338955 ],1e-10);


