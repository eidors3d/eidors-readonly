function [pp, m] = piece_poly_fit(points, fstr, linear_frac)
% PIECE_POLY_FIT: piecewise polynomial fitting toolset
%
% [pp m] = pp_fit(points, fstr) fits a piecewise polynomial to a contour
%    points - [x y] contour to be fitted
%    fstr   - function to use: 'spline' or 'pchip' (default)
%    m      - returns the mean of the original points which was subtracted
%             for fitting
%
% [th xy] = pp_fit(pp, start_th, linear_frac) returns points at linear_frac
% distance along the contour using start_th as the starting angle.
%    pp          - piecewise polynomial structure
%    start_th    - starting angle
%    linear_frac - vector of length fractions (0 to 1) to calculate points
%    th          - angles of the desired points (-pi to pi) (use ppval to
%                  obtain the radii.
%    xy          - cartesian coordinates of the points

% (C) Bartlomiej Grychtol, 2010. Licenced under GPL v2 or v3
% $Id$

if isstr(points) && strcmp(points,'UNIT_TEST'); do_unit_test; return ; end

if nargin < 2
    fhandle = 'pchip';
end
m = [];
if isstruct(points)
    pts = max(100,length(linear_frac)*3);
    start_th = fstr;
    [pp m] = path_len(points, pts, star_th, linear_frac );
else
    [pp m] = fit_to_pp(points, fstr);
end

%%
function [pp centroid] = fit_to_pp(points, fstr)

% 0. Subtract the centroid and convert to polar coords
centroid = mean(points);
n_points = size(points,1);
points = points - repmat(centroid, [n_points,1]);
[ppoints(:,1), ppoints(:,2)] = cart2pol(points(:,1), points(:,2));
% 1. close the loop:
ppoints = sortrows(ppoints,1);
r = ppoints(:,2);rho = ppoints(:,1);
% add points at +/- pi (if none present)
if rho(1) == -pi
    r = [r; r(1)];
    rho = [rho; pi];
elseif rho(end) == pi
    r = [r(end); r];
    rho = [-pi; rho];
else
    dist = 2*pi - rho(end) + rho(1);
    m = r(end)*(pi+rho(1))/dist + r(1)*(pi-rho(end))/dist;
    r = [m; r; m];
    rho = [-pi; ppoints(:,1); pi  ];
end
% 2. fit
switch fstr
    case 'pchip'
        pp=pchip(rho,r);
    case 'spline'
        df = (r(2) - r(end-1)) / ( rho(end-1) + rho(2));
        pp=spline(rho, [df;r;df]);
%         pp=spline(rho, r);
end

function r = fit_from_pp(pp,rho)
r = ppval(pp,rho);

% start_th is starting angle for interpolation
% linear_frac is length fraction at which to find the theta => [0.1, 0.5,]
function [th_frac xy]  = path_len( pp, pts, start_th, linear_frac )
   th = linspace(start_th, start_th+2*pi,pts+1)';
   [x,y] = pol2cart(th,ones(size(th)));
   th = cart2pol(x,y);
   th(end) = [];
   th = sortrows(th);
   dth= diff(th);
   rad = fit_from_pp(pp, th);
   drad= diff(rad);
   rad_= 0.5*(rad(1:end-1) + rad(2:end));
   dlen= sqrt( (rad_ .* dth).^2 + drad.^2 );
   pathlen = [0;cumsum(dlen)];
   
   npath = pathlen/max(pathlen);
   linear_frac = linear_frac + 0.5 + start_th / (2*pi) ;
   idx = find(linear_frac > 1);
   linear_frac(idx) = linear_frac(idx) - 1;
   th_frac = interp1(npath, th, linear_frac);

   lrad = fit_from_pp(pp,th_frac);
   [xe,ye]= pol2cart( th_frac, lrad);
   xy = [xe ye];
%%
function do_unit_test
%     a = [
%    -0.8981   -0.7492   -0.2146    0.3162    0.7935    0.9615    0.6751    0.0565   -0.3635   -0.9745
%     0.1404    0.5146    0.3504    0.5069    0.2702   -0.2339   -0.8677   -0.6997   -0.8563   -0.4668 ]';
% 
% centroid = mean(a);
% n_points = size(a,1);
% a = a - repmat(centroid, [n_points,1]);

th=linspace(0,2*pi,33)'; th(end)=[];
a=[sin(th)*0.3,cos(th)];


  pp= piece_poly_fit(a,'pchip');
   
   th2 = linspace(-pi,pi,33)';th2(end)=[];
   r2 = ppval(pp,th2);
   [xf,yf] = pol2cart(th2,r2);


   
   [lfrac xy] = path_len( pp, 100, 0, linspace(.5,.8,6)' );
   
   clf
   subplot 121
   plot(a(:,1),a(:,2),'*',xf,yf,'b-', xy(:,1),xy(:,2),'r+');
   axis equal

   pp= piece_poly_fit(a,'spline');
   th2 = linspace(-pi,pi,33)';th2(end)=[];
   r2 = ppval(pp,th2);
   [x,y] = pol2cart(th2,r2);
   subplot 122
   plot(a(:,1),a(:,2),'r+'); hold on
   plot(x,y,'.');
   axis equal
