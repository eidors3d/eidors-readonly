% Compare 2D algorithms
% $Id: tutorial120a.m,v 1.1 2007-07-12 18:08:41 aadler Exp $

imb=  mk_common_model('c2c',16);

e= size(imb.fwd_model.elems,1);
bkgnd= 1;

% Solve Homogeneous model
img= eidors_obj('image','');
img.elem_data= bkgnd * ones(e,1);
img.fwd_model= imb.fwd_model;
vh= fwd_solve( img );

% Add Two triangular elements
img.elem_data([25,37,49:50,65:66,81:83,101:103,121:124])=bkgnd * 2;
img.elem_data([95,98:100,79,80,76,63,64,60,48,45,36,33,22])=bkgnd * 2;
vi= fwd_solve( img );

% Add 12dB noise
vi_n= vi; 
vi_n.meas = vi.meas + std(vi.meas - vh.meas)/10^(12/20) ...
                     *randn(size(vi.meas));

show_fem(img);
axis square; axis off
print -r50 -dpng tutorial120a.png;
