function inhomg_img = create_inclusion( homg_img, center, r, inclusion_mat, order)
% USAGE: inhomg_img = create_inclusion( homg_img, inclusion_xyz, ...
%                                       inclusion_r, inclusion_mat, ...
%                                       order)
%
% Parameters: 
%      homg_img  -  image structure (fwd_model + elem_data)
%      center    -  center of spherical inclusion
%      r         -  radius of inclusion
%      inclusion_mat   -  inclusion material
%      order     -  interpolation order (default: 0)
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
   if nargin == 4
     order = 0;
   end
   
   % Calculate interpolation points 
   mdl_pts = interp_mesh(fwd_model, order);
   points_in_element = size(mdl_pts, 3);
  
   % Main loop
   for ei=1:size(mdl_pts,1)
      points_in_sphere = 0;
      for vi=1:points_in_element
         if norm(mdl_pts(ei,:,vi) - center) < r
             points_in_sphere = points_in_sphere + 1;
         end
      end
      ratio = points_in_sphere/points_in_element;
      mat(ei) = ratio*inclusion_mat + (1-ratio)*mat(ei);     
   end

   % Create new image object
   inhomg_img= eidors_obj('image', name, ...
                          'elem_data', mat, ...
                          'fwd_model', fwd_model );
    
  
  
  