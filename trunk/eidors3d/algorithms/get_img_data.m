function [img_data]= get_img_data(img, flag);
% GET_IMG_DATA: get parameter data from eidors image object
% [img_data]= get_img_data(img, flag);
% img_data - img parameter data mapped onto chosen mesh
%
% FLAGS: flag = 0 (default)
%    get data mapped onto elems in the fwd_model

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: get_img_data.m,v 1.1 2008-03-16 11:08:04 aadler Exp $

if nargin==1; flag=0; end

if ~strcmp( img.type, 'image' );
    error('get_img_data expects an image type object');
    return;
end

try
   img_data= img.node_data;
catch
   img_data= img.elem_data;
end

try
   c2f= img.fwd_model.coarse2fine;
catch
   return;
end

if size(img_data,1)==size(c2f,2)
   img_data=c2f * img_data;
end
   
