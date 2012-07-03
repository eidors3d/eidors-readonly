load montreal_data_1995
img= inv_solve(imdl, zc_resp(:,1), zc_resp);
img.calc_colours.ref_level=0;
show_slices(img);
print_convert adult_ex03b.png
img.elem_data= img.elem_data(:,22);
show_fem(img);
print_convert adult_ex03a.png
