function inhomg_img = create_inclusion( homg_img, center, r, inclusion_mat)
% USAGE: inhomg_img = create_inclusion( homg_img, inclusion_xyz, ...
%                                           inclusion_r, inclusion_mat)
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
   dim = size(fwd_model.nodes, 2);
  
   for ei=1:size(fwd_model.elems)
      elem_center = zeros(1,dim);
      for vi=1:(dim+1)
         elem_center = elem_center + ...
         fwd_model.nodes(fwd_model.elems(ei,vi),:);
      end
      elem_center = elem_center/(dim+1);
      if norm(elem_center - center) < r
         mat(ei) = inclusion_mat;     
      end
   end

   inhomg_img= eidors_obj('image', name, ...
                          'elem_data', mat, ...
                          'fwd_model', fwd_model );
    
  end
  
  