% Compare 3D algorithms
% $Id$

imb=  mk_common_model('n3r2',[16 2]);
bkgnd= 1;

% Homogenous Data
img= mk_image(imb.fwd_model, bkgnd);
vh= fwd_solve( img );

% Inhomogenous Data - Load from file 'datacom'
load datacom A B;
img.elem_data(A)= bkgnd*1.2;
img.elem_data(B)= bkgnd*0.8;
clear A B;
vi= fwd_solve( img );

% Add 15dB noise
vi_n= vi; 
vi_n.meas = vi.meas + std(vi.meas - vh.meas)/10^(10/20) ...
                     *randn(size(vi.meas));
sig= sqrt(norm(vi.meas - vh.meas));

subplot(121);
show_fem(img); axis square;

subplot(122);
show_fem(img); axis square;
crop_model([],  inline('y<0','x','y','z'))
view(-51,14);
print_convert('tutorial130a.png', '-density 100')
