% fmdl= ng_mk_cyl_models(3,[7,1],[0.2,0,0.05]); 
pl= 1; for SKIP = [1,2,6];
    subplot(1,3,pl); pl=pl+1;

imdl = mk_common_model('f2c2',16);
fmdl = imdl.fwd_model;
fmdl.stimulation = mk_stim_patterns(16,1,[0,SKIP],[0,1],{},1);

img = mk_image(fmdl,1);
J = calc_jacobian(img)';
sens = img;
s1 = zeros(size(J,1),1);
for i=1:size(J,2)
    s1 = s1 + J(:,i).^2;
end
sens.elem_data = sqrt(s1)./get_elem_volume(fmdl);

sens.calc_colours.clim = 3;
sens.calc_colours.ref_level = 0;
%show_fem(sens,[0 1]);
show_slices(sens)

end
%eidors_colourbar(sens);