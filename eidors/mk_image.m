function img= mk_image( elem_data, varargin)
% MK_IMAGE: create eidors_image object
% Utility function to create an eidors_image object:
%   img= mk_image( elem_data, param1, value1, ...)
% Usage: Default, create image with just elem_data
%  img = mk_image( elem_data );
% 
% Usage: create image with other field values
%  img = mk_image( elem_data, 'name', 'image name', 'param2', 'value' );
%
% Usage: create image from previous image, override fields
%  img = mk_image( other_image, 'name', 'image name', 'param2', 'value' );

% (C) 2008 Andy Adler. Licenced under GPL version 2 or 3
% $Id$

if isstruct(elem_data)
   img= elem_data;
else
   img= eidors_obj('image','unnamed','elem_data',elem_data);
end

for i=1:2:length(varargin)
   img.( varargin{i} ) = varargin{i+1};
end
