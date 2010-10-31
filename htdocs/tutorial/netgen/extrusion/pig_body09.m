imdl = mk_common_model('c2c2',16);
imdl.fwd_model.electrode = imdl.fwd_model.electrode([1,16:-1:2]);
imdl.fwd_model.normalize_measurements= 1;
imr= inv_solve(imdl, vh, vi);

show_fem(imr); axis tight; axis off

print_convert pig_body09a.jpg

