% fwd_model $Id$

fmdl = mk_library_model('cylinder_16x1el_fine');

pixel_grid= 32;
nodes= fmdl.nodes;
xyzmin= min(nodes,[],1);  xyzmax= max(nodes,[],1);
xvec= linspace( xyzmin(1), xyzmax(1), pixel_grid+1);
yvec= linspace( xyzmin(2), xyzmax(2), pixel_grid+1);
zvec= [0.6*xyzmin(3)+0.4*xyzmax(3), 0.4*xyzmin(3)+0.6*xyzmax(3)];

% CALCULATE MODEL CORRESPONDENCES
[rmdl,c2f] = mk_grid_model(fmdl, xvec, yvec, zvec);


% SHOW MODEL CORRESPONDENCE

clf;
show_fem(fmdl);  % fine model
crop_model(gca, inline('x-z<-8','x','y','z'))
hold on
hh=show_fem(rmdl); set(hh,'EdgeColor',[0,0,1]);
hold off

view(-47,18); 

print_convert framework01a.png
