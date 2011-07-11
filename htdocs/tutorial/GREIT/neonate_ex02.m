vh = mean(vv,2);        % reference is average
vi = vv(:,[45,70,173]); %3 inspirations

img = inv_solve(imdl,vh,vi);

img.show_slices.img_cols = 3;
img.show_slices.sep      = 2;
img.calc_colours.ref_level=0;
show_slices(img);

print_convert neonate_ex02a.png
