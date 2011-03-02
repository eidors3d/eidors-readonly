function img= mk_image(mdl, elem_data, name)
% MK_IMAGE: create eidors image object
%   img= mk_image(mdl, elem_data, name)
%
% Utility function to create an eidors_image object:
% Usage 1:
%   img = mk_image( inv_model )  -> uses jacobian backgnd for conductivity
%
% Usage 2: mdl can be a fwd_model or inv_model
%   img = mk_image( mdl, 1 ) -> uniform image with conductivity 1
%   img = mk_image( mdl, 2*ones(n_elems,1) ) -> uniform with c=2
%   img = mk_image( mdl, 1, 'This name') -> Specify a 'name' attribute
% 
% Usage 3: create image from previous image, override conductity
%  img = mk_image( other_image, 2 ) -> image with c=2

% (C) 2008-10 Andy Adler. Licenced under GPL version 2 or 3
% $Id$

if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end

if nargin<3; name = 'Created by mk_image'; end
if nargin<2; 
   try
     elem_data = mdl.jacobian_bkgnd.value;
   catch
     error('mk_image: for one parameter, needs a mdl.jacobian_bkgnd field');
   end
end

switch mdl.type
   case 'inv_model'; mdl = mdl.fwd_model;
   case 'fwd_model'; mdl = mdl; % keep model
   case 'image';     mdl = mdl.fwd_model;
   otherwise; error('mk_image: no inv_model, fwd_model or image object');
end

img = eidors_obj('image',name);
img.fwd_model = mdl;
img.elem_data = NaN*ones(size(mdl.elems,1),1);
img.elem_data(:) = elem_data;

% TESTS:
function do_unit_test
   imdl = mk_common_model('a2c2',8);
   im0 = mk_image( imdl );
   unit_test_cmp('im1',im0.elem_data, ones(64,1) );

   im0 = mk_image( imdl, 2 );
   unit_test_cmp('im2',im0.elem_data, 2*ones(64,1) );

   im0 = mk_image( imdl.fwd_model, 3*ones(64,1) );
   unit_test_cmp('im3',im0.elem_data, 3*ones(64,1) );
