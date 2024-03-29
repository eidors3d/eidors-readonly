% Simulate Moving Ball - Helix $Id$

% get ng_mdl_16x2_vfine from data_contrib section of web page
n_sims= 20;
stim = mk_stim_patterns(16,2,'{ad}','{ad}',{},1);
fmdl = mk_library_model('cylinder_16x2el_vfine');
fmdl.stimulation = stim;
[vh,vi,xyzr_pt]= simulate_3d_movement( n_sims, fmdl);

clf; show_fem(fmdl)
crop_model(gca, inline('x-z<-15','x','y','z'))
view(-23,14)

hold on
[xs,ys,zs]=sphere(10);
for i=1:n_sims
   xp=xyzr_pt(1,i); yp=xyzr_pt(2,i);
   zp=xyzr_pt(3,i); rp=xyzr_pt(4,i);
   hh=surf(rp*xs+xp, rp*ys+yp, rp*zs+zp);
   set(hh,'EdgeColor',[.4,0,.4],'FaceColor',[.2,0,.2]);
end
hold off

print_convert centre_slice01a.png '-density 100'
