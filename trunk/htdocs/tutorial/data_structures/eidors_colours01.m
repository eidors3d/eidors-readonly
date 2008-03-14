% Show EIDORS colours $Id: eidors_colours01.m,v 1.1 2008-03-14 15:26:02 aadler Exp $

% Create sample image
imdl= mk_common_model('a2c2',16);
img= eidors_obj('image','small demo','fwd_model',imdl.fwd_model, ...
                'elem_data', zeros(1,64) );
img.elem_data(1:2)=[1,1];
img.elem_data([14,16])=[-1,-1];
img.elem_data([27,30])=.5*[1,1];
idx= [51:52, 56:-1:53, 57:60, 64:-1:61, 49,50]; % clockwise index elements 
img.elem_data(idx)= linspace(-2,2,16);

img.elem_data = img.elem_data + 10;
