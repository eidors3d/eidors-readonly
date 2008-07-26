% fwd_model $Id$
load ng_mdl_16x1_fine

pixel_grid= 32;
nodes= ng_mdl_16x1_fine.nodes;
xyzmin= min(nodes,[],1);  xyzmax= max(nodes,[],1);
xvec= linspace( xyzmin(1), xyzmax(1), pixel_grid+1);
yvec= linspace( xyzmin(2), xyzmax(2), pixel_grid+1);
zvec= [0.6*xyzmin(3)+0.4*xyzmax(3), 0.4*xyzmin(3)+0.6*xyzmax(3)];

% CALCULATE MODEL CORRESPONDENCES
[rmdl,c2f] = mk_grid_model(ng_mdl_16x1_fine, xvec, yvec, zvec);


% SHOW MODEL CORRESPONDENCE

clf;
show_fem(ng_mdl_16x1_fine);  % fine model
crop_model(gca, inline('x-z<-8','x','y','z'))
hold on
show_fem(rmdl);
hold off

view(-47,18); 

print -dpng -r100 framework01a.png
