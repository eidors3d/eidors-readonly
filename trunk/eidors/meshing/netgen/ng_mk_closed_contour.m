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

[ellip_shape,params] = proc_shape(ellip_shape);
[fmdl, mat_idx] = ng_mk_ellip_models(ellip_shape, elec_pos, ...
                  elec_shape, extra_ng_code);
fmdl.nodes(:,1:2) = fmdl.nodes(:,1:2)*params.rot;

function [ellip_shape,params] = proc_shape(ellip_shape);
  height = ellip_shape{1}; 
  xy_shape=ellip_shape{2};
  if length(ellip_shape)>=3
    maxh= ellip_shape{3};
  else
    maxh= [];
  end

  [u,d,v]= svd(cov(xy_shape));
  params.rot = [1,0;0,-1]*u;
  ellip_ax = sqrt(2*[d(1,1),d(2,2)])  
pause(1);
  ellip_shape= [height, ellip_ax, maxh];
