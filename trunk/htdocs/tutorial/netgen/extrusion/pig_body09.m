imdl = mk_common_model('c2c2',16);
imdl.fwd_model.electrode = imdl.fwd_model.electrode([1,16:-1:2]);
imdl.fwd_model = mdl_normalize(imdl.fwd_model, 1);
imr= inv_solve(imdl, vh, vi);
clf
show_fem(imr); axis tight; axis off; axis equal

print_convert pig_body09a.jpg

