function inhomg_img = create_inclusion( homg_img, center, r, inclusion_mat, order)
% USAGE: inhomg_img = create_inclusion( homg_img, inclusion_xyz, ...
%                                       inclusion_r, inclusion_mat )
%
% Parameters: 
%      homg_img  -  image structure (fwd_model + elem_data)
%      center    -  center of spherical inclusion
%      r         -  radius of inclusion
%      inclusion_mat   -  inclusion material
%
% Function create inclusion by changing material properties for elements
% which center is is inside of the sphere defined by (center and r).
%
% (C) 2009,  Bartosz Sawicki
% $Id$
% Licenced under the GPLv2 or later

name = homg_img.name;
mat = homg_img.elem_data;
fwd_model = homg_img.fwd_model;

% Use coarse2fine mapping to determine if elements are in sphere
xyzr = [center r];
fraction = mk_c2f_circ_mapping( homg_img.fwd_model, xyzr');

% Inclusion is overlaped on old material
mat = fraction*inclusion_mat + (1-fraction).*mat;     

% Create new image object
inhomg_img= eidors_obj('image', name, ...
                       'elem_data', mat, ...
                       'fwd_model', fwd_model );
    
  
  
  