% Create 3D FEM model of the gallery
% $Id: tutorial410a.m,v 1.1 2007-06-15 18:59:35 aadler Exp $
n_rings= 9;
factor= 2;
levels= [-6 -4 -2.5 -1.5 -1 -0.5 -0.25 0 0.25 0.5 1 1.5 2.5 4 6];

Electrode_Positions_Ring1_EZG04;
elec_posn= EZG04_Ring1;

Anneau1_Juillet2004_wen32_1;
data_tomel= Data_Ring1_July2004_Wen32_1;
real_data= mk_data_tomel(data_tomel,'Mont-Terri data','Wenner protocol');

gallery_3D_fwd = mk_gallery(elec_posn,data_tomel,n_rings,factor,levels);

subplot(121)
show_fem(gallery_3D_fwd); axis square; view(0.,15.);
subplot(122)
show_fem(gallery_3D_fwd); axis square; view(0.,75);
%print -r100 -dpng tutorial410a.png;

