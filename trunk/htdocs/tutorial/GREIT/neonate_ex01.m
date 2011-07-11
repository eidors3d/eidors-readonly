imdl = mk_common_gridmdl('GREITc1');

% Data: eidors3d.sf.net/data_contrib/if-neonate-spontaneous/if-neonate-spontaneous.zip
vv= eidors_readdata('P04P-1016.get');
vh = mean(vv,2);
vi = vv(:,[45,70,173]); %3 inspirations
img = inv_solve(imdl,vh,vi);
img.show_slices.img_cols = 3;
img.show_slices.sep      = 2;
img.calc_colours.ref_level=0;
show_slices(img);

return
img = inv_solve(imdl,vh,vv);
plot(sum(img.elem_data));
