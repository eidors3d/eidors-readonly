% TV: Create object $Id: TV_hyperparams01.m,v 1.1 2008-03-14 15:57:31 aadler Exp $

imb=  mk_common_model('c2c2',16);
fmdl= imb.fwd_model;
sigma= ones(size(fmdl.elems,1),1);
img= eidors_obj('image','','elem_data',sigma,'fwd_model',fmdl);

vh= fwd_solve( img );

sigma([25,37,49:50,65:66,81:83,101:103,121:124])=2;
sigma([95,98:100,79,80,76,63,64,60,48,45,36,33,22])=2;
    
img.elem_data= sigma;
vi= fwd_solve( img );

sig= sqrt(norm(vi.meas - vh.meas));
m= size(vi.meas,1);
vi.meas = vi.meas + .0001*sig*randn(m,1);

subplot(221)
show_fem(img);
axis equal; axis off;
print -r100 -dpng TV_hyperparams01a.png;
