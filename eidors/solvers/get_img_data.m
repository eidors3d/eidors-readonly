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
if ischar(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

if ~strcmp( img(1).type, 'image' );
    error('get_img_data expects an image type object');
    return;
end

if ~isfield(img,'elem_data')
   img = data_mapper(img);
end

if isfield(img(1),'elem_data')
   for i=1:length(img(:))
      try 
         img_data{i}= [img(i).elem_data(:,img(i).get_img_data.frame_select)];
      catch
         % concatenate img.elem_data for each img
         img_data{i}= [img(i).elem_data];
      end
   end
   img_data = [img_data{:}];
elseif isfield(img,'node_data')
   for i=1:length(img(:))
      try 
         img_data{i}= [img(i).node_data(:,img(i).get_img_data.frame_select)];
      catch
         img_data{i}= [img(i).node_data];
      end
   end
   img_data = [img_data{:}];
else
   error('No field elem_data or node_data found on img');
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
   unit_test_cmp('elem_02', get_img_data(imgk), ones(64,1) )


   imgk = img; imgk(2) = img;
   unit_test_cmp('elem_03', get_img_data(imgk), ones(64,2) )

   % EXPECTED FAIL
   % imgk(2).elem_data = imgk(2).elem_data';

   imgk(1).elem_data = imgk(1).elem_data*[1,2];
   unit_test_cmp('elem_04', get_img_data(imgk(1)), ones(64,1)*[1,2] )
   unit_test_cmp('elem_05', get_img_data(imgk), ones(64,1)*[1,2,1] )

   imgk = img; imgk.elem_data = img.elem_data*[1,2];
   imgk.get_img_data.frame_select = 2;
   unit_test_cmp('elem_06', get_img_data(imgk), 2*ones(64,1) )
   % TODO: This fails because we need a loop in the test function.
   %   Lookup whether this can be vectorized in matlab

   imgk(2) = imgk(1); imgk(2).elem_data = img.elem_data*[3,4];
   imgk(2).get_img_data.frame_select = 2;
   unit_test_cmp('elem_07', get_img_data(imgk), ones(64,1)*[2,4] )
   imgk(2).get_img_data.frame_select = 1:2;
   unit_test_cmp('elem_08', get_img_data(imgk), ones(64,1)*[2,3,4] )

   nd =ones(num_nodes(img),1);
   imgn = rmfield(img,'elem_data'); imgn.node_data=nd;
   unit_test_cmp('node_01', get_img_data(imgn), nd)

   imgn.node_data = nd*[3,4,5];
   unit_test_cmp('node_02', get_img_data(imgn), nd*[3,4,5]);
   imgn(2) = imgn;
   unit_test_cmp('node_03', get_img_data(imgn), nd*[3,4,5,3,4,5]);
   imgn(1).get_img_data.frame_select = 1;
   unit_test_cmp('node_04', get_img_data(imgn), nd*[3,3,4,5]);
   imgn(2).get_img_data.frame_select = 3;
   unit_test_cmp('node_05', get_img_data(imgn), nd*[3,5]);
   

   i2 = mk_common_model('b2c0',8); fmdl2= i2.fwd_model;
   c2f = mk_coarse_fine_mapping(fmdl2, img.fwd_model);
