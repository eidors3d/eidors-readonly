load montreal_data_1995
img= inv_solve(imdl, zc_resp(:,1), zc_resp);

img.show_slices.img_cols=9;
show_slices(img);
print_convert adult_ex03b.png

img.calc_colours.ref_level = 0;
img.get_img_data.frame_select = 22;
show_fem(img);
print_convert adult_ex03a.png
