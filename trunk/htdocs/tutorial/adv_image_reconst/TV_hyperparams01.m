% TV: Create object $Id$

imb=  mk_common_model('c2c2',16);
img= mk_image(imb.fwd_model, 1);

vh= fwd_solve( img );

img.elem_data([25,37,49:50,65:66,81:83,101:103,121:124])=2;
img.elem_data([95,98:100,79,80,76,63,64,60,48,45,36,33,22])=2;
    
vi= fwd_solve( img );

subplot(221)
show_fem(img);
axis equal; axis off;
print_convert TV_hyperparams01a.png '-density 100';
