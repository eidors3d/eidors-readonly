% fwd_model $Id$

fmdl = mk_library_model('cylinder_16x1el_fine');
fmdl.stimulation = mk_stim_patterns(16,1,'{ad}','{ad}',{},1);

subplot(121)
show_fem(fmdl);
view([-5 28]);
crop_model(gca, inline('x+2*z>20','x','y','z'))

subplot(122)
show_fem(fmdl);
view(0,0);
crop_model(gca, inline('y>-10','x','y','z'))
set(gca,'Xlim',[-4,4],'Zlim',[-2,2]+5);

print_convert framework00a.png
