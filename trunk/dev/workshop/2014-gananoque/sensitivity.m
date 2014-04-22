% fmdl= ng_mk_cyl_models(3,[7,1],[0.2,0,0.05]); 
imdl = mk_common_model('f2c2',16);
fmdl = imdl.fwd_model;
show_fem(fmdl);
img = mk_image(fmdl,1);
J = calc_jacobian(img)';
sens = img;
sens.elem_data = J(:,5)./get_elem_volume(fmdl);

show_fem(sens);
eidors_colourbar(sens);