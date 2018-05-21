fname = 'DATA/STUDYNAME/SUBJECT_1/YYYYMMDD/Eit/Viasys/1001_c4.get';
vv= eidors_readdata(fname);
img= inv_solve(imdl, mean(vv,2), vv);

imgs= -calc_slices(img); % Negative to air is +
imgs(isnan(imgs(:)))= 0;

img.calc_colours.ref_level=0;
img.elem_data = img.elem_data(:,2:4:120);
img.show_slices.img_cols = 10;
clf; show_slices(img);

print_convert 'GREIT_IBEX_03a.jpg'
