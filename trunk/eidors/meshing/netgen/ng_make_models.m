function fmdl = ng_make_cyl_models(cyl_shape, elec_pos, elec_shape, extra_geo_code);
% NG_MAKE_CYL_MODELS: create cylindrical models using netgen
% 
% cyl_shape = {height, [radius, [minsz]]}
%    if height = 0 -> calculate a 2D shape
%       radius     -> default = 1
%       minsz      -> coarse mesh (netgen default)
%
% elec_pos = [n_elecs_per_plane,z_planes] 
%    OR
% elec_pos = [x,y,z] centres

%{[width,height,minsz]}}
% Usage:
%  fmdl= ng_make_cyl_models( 

% (C) Andy Adler, 2009. Licenced under GPL v2 or v3
% $Id$

