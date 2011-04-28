load montreal_data_1995
img= inv_solve(imdl, zc_resp(:,22), zc_resp);
show_fem(img);
print_convert adult_ex03a.png
show_slices(img);
print_convert adult_ex03b.png
