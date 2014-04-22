% fmdl= ng_mk_cyl_models(3,[7,1],[0.2,0,0.05]); 
slices= [];
for pl=1:3
    subplot(2,2,pl); 
    switch pl;
        case 1; CSKIP = 1; MSKIP= 1;
        case 2; CSKIP = 6; MSKIP= 1;
        case 3; CSKIP = 6; MSKIP= 6;
    end

imdl = mk_common_model('f2c2',16);
fmdl = imdl.fwd_model;
fmdl.stimulation = mk_stim_patterns(16,1,[0,CSKIP],[0,MSKIP],{},1);

img = mk_image(fmdl,1);
J = calc_jacobian(img)';
sens = img;
s1 = zeros(size(J,1),1);
for i=1:size(J,2)
    s1 = s1 + J(:,i).^2;
end
sens.elem_data = sqrt(s1)./get_elem_volume(fmdl);

%sens.calc_colours.clim = 3;
sens.calc_colours.ref_level = 0;
%show_fem(sens,[0 1]);
show_slices(sens)
eidors_colourbar(sens);
sli = calc_slices(sens);
slices = [slices,sli(:,32)./max(sli(:,32))];
end
%eidors_colourbar(sens);
subplot(2,2,4)
plot(slices);