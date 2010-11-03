% Show images $Id$

subplot(121)
axis square; view(30.,80.);
show_fem(gallery_3D_img);

subplot(122)
gallery_3D_resist= gallery_3D_img; % Create resistivity image
gallery_3D_resist.elem_data= 1./gallery_3D_img.elem_data;
show_slices(gallery_3D_resist,[inf,inf,0]);

%print -r125 -dpng tutorial410d.png;
