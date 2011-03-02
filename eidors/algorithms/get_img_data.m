function [img_data, n_images]= get_img_data(img)
% GET_IMG_DATA: get parameter data from eidors image object
% [img_data, n_images]= get_img_data(img)
% img_data - img parameter data mapped onto chosen mesh
%   get_img_data looks at the elem_data or node_data
%   parameters in the image. If a course2fine parameter
%   exists in img.fwd_model, then it is used. 

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

% TODO
% FLAGS: flag = 0 (default)
%    get data mapped onto elems in the fwd_model
%if nargin==1; flag=0; end

% TEST
if isstr(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

if ~strcmp( img(1).type, 'image' );
    error('get_img_data expects an image type object');
    return;
end

try
   img_data= [img.elem_data];
catch
   img_data= [img.node_data];
end

if size(img_data,1)==1; img_data=img_data';end

n_images= size(img_data,2);

try
   c2f= img.fwd_model.coarse2fine;
catch
   return;
end

if size(img_data,1)==size(c2f,2)
   img_data=c2f * img_data;
end

% TESTS:
function do_unit_test
   img=calc_jacobian_bkgnd(mk_common_model('a2c0',8)); 
   unit_test_cmp('elem_01', get_img_data(img), ones(64,1) )

   imgk = img; imgk.elem_data = img.elem_data';
   unit_test_cmp('elem_02', get_img_data(imgk), ones(64,2) )

   imgk = img; imgk(2) = img;
   unit_test_cmp('elem_03', get_img_data(imgk), ones(64,2) )

   % EXPECTED FAIL
   % imgk(2).elem_data = imgk(2).elem_data';

   imgk(1).elem_data = imgk(1).elem_data*[1,2];
   unit_test_cmp('elem_04', get_img_data(imgk(1)), ones(64,1)*[1,2] )
   unit_test_cmp('elem_05', get_img_data(imgk), ones(64,1)*[1,2,1] )

   nd =ones(size(img.fwd_model.nodes,1),1);
   imgn = rmfield(img,'elem_data'); imgn.node_data=nd;
   unit_test_cmp('node_01', get_img_data(imgn), nd)

   i2 = mk_common_model('b2c0',8); fmdl2= i2.fwd_model;
   c2f = mk_coarse_fine_mapping(fmdl2, img.fwd_model);
   img.fwd_model.coarse2fine = c2f;
   unit_test_cmp('c2f_01', get_img_data(img), c2f*ones(64,1))
    
   % TODO ADD TESTS WHICH GIVE C2F ON NODE_DATA
