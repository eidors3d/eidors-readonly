% Four cut planes
cut_planes = -[0.6;0.7;0.8;0.9];
imgv.calc_colours.backgnd = 0.9*[1,1,1];
show_slices(imgv,cut_planes*[inf,1,inf]);
print_convert view_3D_surf03a.png '-density 150'
