function [img_data, n_images]= get_img_data(img, flag);
% GET_IMG_DATA: get parameter data from eidors image object
% [img_data, n_images]= get_img_data(img, flag);
% img_data - img parameter data mapped onto chosen mesh
% n_images - number of images in img
%
% FLAGS: flag = 0 (default)
%    get data mapped onto elems in the fwd_model

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: get_img_data.m,v 1.3 2008-03-30 21:48:34 aadler Exp $

if nargin==1; flag=0; end

if ~strcmp( img(1).type, 'image' );
    error('get_img_data expects an image type object');
    return;
end

try
   img_data= [img.node_data];
catch
   img_data= [img.elem_data];
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

