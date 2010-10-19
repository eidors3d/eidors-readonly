function C = fourier_fit(points,N);
% FOURIER_FIT: fit fourier descriptors to a contour

% (C) Andy Adler, 2010. Licenced under GPL v2 or v3
% $Id$

% Fourier Fit
% [1, cos(t1), cos(2*t1), ... sin(t1) ...] * [C1] = [R1]
%            ...                         
% [1, cos(tN), cos(2*tN), ... sin(tN) ...] * [CN] = [RM]

if isstr(points) && strcmp(points,'UNIT_TEST'); do_unit_test; return ; end

C = fit_to_fourier(points,N);
    
function C = fit_to_fourier(xy_shape,N);
   xc = xy_shape(:,1);
   yc = xy_shape(:,2);

   ang = atan2(yc,xc);
   rad = sqrt(yc.^2 + xc.^2);
   A= fit_matrix(ang, N);
   AR = diag(diag(A'*A));
   C = (A'*A + 1e-9*AR)\A'*rad;
   
function R = fit_from_fourier(ang,C);
   A = fit_matrix(ang, (length(C)-1)/2);
   R= A*C;

function A= fit_matrix(ang, N);
   A = ones(length(ang), 2*N+1);
   for i=1:N
     A(:,i+1  ) = cos(i*ang);
     A(:,i+1+N) = sin(i*ang);
   end

% start_th is starting angle for interpolation
% linear_frac is length fraction at which to find the theta => [0.1, 0.5,]
function [pathlen, th_frac]  = path_len( C, pts, start_th, linear_frac )
   th = linspace(start_th, start_th+2*pi,pts+1)';
   dth= diff(th);
   rad = fit_from_fourier(th, C);
   drad= diff(rad);
   rad_= 0.5*(rad(1:end-1) + rad(2:end));
   dlen= sqrt( (rad_ .* dth).^2 + drad.^2 );
   pathlen = [0;cumsum(dlen)];
   
   npath = pathlen/max(pathlen);
   th_frac = interp1(npath, th, linear_frac);
   
function do_unit_test

    a = [
   -0.8981   -0.7492   -0.2146    0.3162    0.7935    0.9615    0.6751    0.0565   -0.3635   -0.9745
    0.1404    0.5146    0.3504    0.5069    0.2702   -0.2339   -0.8677   -0.6997   -0.8563   -0.4668 ]';

th=linspace(0,2*pi,33)'; th(end)=[];
a=[sin(th)*0.3,cos(th)];


   C= fourier_fit(a,6);
   ang = linspace(0,2*pi,50)';
   rad = fit_from_fourier(ang,C);

   [xf,yf]= pol2cart( ang, rad);

   [pathlen, lfrac] = path_len( C, 100, 0, linspace(.5,.8,6)' );
   lrad = fit_from_fourier(lfrac,C);
   [xe,ye]= pol2cart( lfrac, lrad);

   plot(a(:,1),a(:,2),'*',xf,yf,'b-', xe,ye,'r+');

