% Compare 2D algorithms
% $Id$

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

% Add -12dB SNR
vi_n= vi; 
nampl= std(vi.meas - vh.meas)*10^(-18/20);
vi_n.meas = vi.meas + nampl *randn(size(vi.meas));

subplot(221)
show_fem(img);
axis square; axis off
print -r100 -dpng tutorial120a.png;
