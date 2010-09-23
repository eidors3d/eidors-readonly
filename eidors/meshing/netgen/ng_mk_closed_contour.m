function [fmdl,mat_idx] = ng_mk_ellip_models(ellip_shape, elec_pos, ...
                  elec_shape, extra_ng_code);
% NG_MK_CLOSED_CONTOUR: fit elliptical model to a contour
%[fmdl,mat_idx] = ng_mk_closed_contour(shape, elec_pos, ...
%                 elec_shape, extra_ng_code);
%
% This function creates netgen models of a closed contour in x,y.
%   An elliptical FEM is created and then perturbed using Fourier 
%   descriptors into the final model.
% It should work well for cases where the contour is roughtly elliptical
%
% INPUT:
% ellip_shape = cell{height, xy_points, [maxsz]]}
%    if height = 0 -> calculate a 2D shape
%    xy_points     -> Npoints x 2. points on the object boundary
%    maxsz  (OPT)  -> max size of mesh elems (default = course mesh)
%
% ELECTRODE POSITIONS:
%  elec_pos = [n_elecs_per_plane,z_planes] 
%     OR
%  elec_pos = [degrees,z] centres of each electrode (N_elecs x 2)
% Note that electrode positions are in the ellipse before fitting
%
% ELECTRODE SHAPES::
%  elec_shape = [width,height, maxsz]  % Rectangular elecs
%     OR
%  elec_shape = [radius, 0, maxsz ]    % Circular elecs
%     OR 
%  elec_shape = [0, 0, maxsz ]         % Point electrodes
%    (point elecs does some tricks with netgen, so the elecs aren't exactly where you ask)
%
% Specify either a common electrode shape or for each electrode
%
% EXTRA_NG_CODE
%   string of extra code to put into netgen geo file. Normally this
%   would be to insert extra materials into the space
%
% OUTPUT:
%  fmdl    - fwd_model object
%  mat_idx - indices of materials (if extra_ng_code is used)
%    Note mat_idx does not work in 2D. Netgen does not provide it.
%
%
% USAGE EXAMPLES:
% Simple 3D ellipse. Major, minor axes = [1.5, 0.8]. No electrodes
%     fmdl= ng_mk_ellip_models([1,1.5,0.8],[0],[]);  show_fem(fmdl);
% 

% (C) Andy Adler, 2010. Licenced under GPL v2 or v3
% $Id$

if nargin < 4; extra_ng_code = {'',''}; end
% Caching should be done in the ng_mk_ellip_models

[ellip_shape,params] = analyse_shape(ellip_shape);
[fmdl, mat_idx] = ng_mk_ellip_models(ellip_shape, elec_pos, ...
                  elec_shape, extra_ng_code);

% Rotate and move the objects
fmdl.nodes(:,1:2) = fmdl.nodes(:,1:2)*params.rot;
fmdl.nodes(:,1) = fmdl.nodes(:,1) + params.xy_mean(1);
fmdl.nodes(:,2) = fmdl.nodes(:,2) + params.xy_mean(2);

function [ellip_shape,params] = analyse_shape(ellip_shape);
  height = ellip_shape{1}; 
  xy_shape=ellip_shape{2};
  if length(ellip_shape)>=3
    maxh= ellip_shape{3};
  else
    maxh= [];
  end

  [u,d,v]= svd(cov(xy_shape));
  params.rot = [1,0;0,-1]*u;
  params.xy_mean = mean(xy_shape,1);

  ellip_ax = sqrt(2*[d(1,1),d(2,2)])  
  ellip_shape= [height, ellip_ax, maxh];

  N = ceil( size(xy_shape,1) / 3 ); % can't fit too many components
  paramc.fourier_fit = fourier_fit(xy_shape, N);

function C = fourier_fit(xy_shape,N);
% Fourier Fit
% [1, cos(t1), cos(2*t1), ... sin(t1) ...] * [C1] = [R1]
%            ...                         
% [1, cos(tN), cos(2*tN), ... sin(tN) ...] * [CN] = [RM]

   xc = xy_shape(:,1); xm= mean(xc);   xc= xc-xm;
   yc = xy_shape(:,2); ym= mean(yc);   yc= yc-ym;

   ang = atan2(yc,xc);
   rad = sqrt(yc.^2 + xc.^2);
   A = ones(length(ang), 2*N+1);
   for i=1:N
     A(:,i+1  ) = cos(i*ang);
     A(:,i+1+N) = sin(i*ang);
   end
   C = A\rad;

function ff2

   ang = linspace(0,2*pi,50)';
   A = ones(length(ang), 2*N+1);
   for i=1:N
     A(:,i+1  ) = cos(i*ang);
     A(:,i+1+N) = sin(i*ang);
   end
   rad = A*C;
   xf = rad.*cos(ang) + xm;
   yf = rad.*sin(ang) + ym;

   if doplot
      hold on;plot(xf,yf,'r-*'); hold off;
   end
