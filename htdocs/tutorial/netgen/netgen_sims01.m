% Netgen simulation $Id$

extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05);'};
fmdl= ng_mk_cyl_models(0,[9],[0.2,0,0.05],extra); 
img= eidors_obj('image','ball'); img.fwd_model= fmdl;
ctr = interp_mesh(fmdl); ctr=(ctr(:,1)-0.2).^2 + (ctr(:,2)-0.2).^2;
img.elem_data = 1 + 0.1*(ctr < 0.2^2);

subplot(221)
show_fem(img);

print -dpng -r100 netgen_sims01a.png
