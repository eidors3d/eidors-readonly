% Netgen simulation $Id$

extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05);'};
fmdl= ng_mk_cyl_models(0,[9],[0.2,0,0.05],extra); 
ctr = interp_mesh(fmdl); ctr=(ctr(:,1)-0.2).^2 + (ctr(:,2)-0.2).^2;
img= mk_image(fmdl, 1 + 0.1*(ctr < 0.2^2));

subplot(221)
show_fem(img);

print_convert netgen_sims01a.png '-density 100'
