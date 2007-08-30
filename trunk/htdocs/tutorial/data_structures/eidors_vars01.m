% Show EIDORS colours $Id: eidors_vars01.m,v 1.3 2007-08-30 03:58:28 aadler Exp $

imdl= mk_common_model('a2c2',16);
img= eidors_obj('image','small demo','fwd_model',imdl.fwd_model, ...
                'elem_data', zeros(1,64) );
img.elem_data(1:2)=[1,1];
img.elem_data([14,16])=[-1,-1];
img.elem_data([27,30])=.5*[1,1];
show_fem(img);

print -r50 -dpng eidors_vars01a.png
