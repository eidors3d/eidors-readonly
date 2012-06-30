% Show EIDORS colours $Id$

% Create sample image
imdl= mk_common_model('a2c2',16);
img= mk_image(imdl.fwd_model, 0);
img.elem_data(1:2)=[1,1];
img.elem_data([14,16])=[-1,-1];
img.elem_data([27,30])=.5*[1,1];
idx= [51:52, 56:-1:53, 57:60, 64:-1:61, 49,50]; % clockwise index elements 
img.elem_data(idx)= linspace(-2,2,16);

img.elem_data = img.elem_data + 10;
