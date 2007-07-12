% Compare 3D algorithms
% $Id: tutorial130a.m,v 1.1 2007-07-12 18:08:41 aadler Exp $

imb=  mk_common_model('n3r2',16);
e= size(imb.fwd_model.elems,1);
bkgnd= 1;

% Homogenous Data
img= eidors_obj('image','Demo Image');
img.elem_data= bkgnd*ones(e,1);
img.fwd_model= imb.fwd_model;
vh= fwd_solve( img );

% Inhomogenous Data - Load from file 'datacom'
load datacom A B;
img.elem_data(A)= bkgnd*1.2;
img.elem_data(B)= bkgnd*0.8;
clear A B;
vi= fwd_solve( img );

% Add 12dB noise
vi_n= vi; 
vi_n.meas = vi.meas + std(vi.meas - vh.meas)/10^(12/20) ...
                     *randn(size(vi.meas));
sig= sqrt(norm(vi.meas - vh.meas));

subplot(121);
show_fem(img); axis square;

subplot(122);
show_fem(img); axis square;
crop_model([],  inline('y<0','x','y','z'))
view(-51,14);
print -r100 -dpng tutorial130a.png;
