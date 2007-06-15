% Show images $Id: tutorial410d.m,v 1.1 2007-06-15 19:56:56 aadler Exp $

subplot(121)
axis square; view(30.,80.);
show_fem(gallery_3D_img);

subplot(122)
gallery_3D_resist= gallery_3D_img; % Create resistivity image
gallery_3D_resist.elem_data= 1./gallery_3D_img.elem_data;
show_slices(gallery_3D_resist,[inf,inf,0],20,20);

 print -r100 -dpng tutorial410d.png;
