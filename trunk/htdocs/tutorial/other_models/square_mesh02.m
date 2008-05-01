% Create and show square model $Id: square_mesh02.m,v 1.3 2008-05-01 16:27:44 aadler Exp $

% Create square mesh model
cmdl= mk_grid_model( linspace(-8,8,17), linspace(-11.5,-0.5,13) );

clf;
show_fem(fmdl);
hold on;
h= trimesh(cmdl.elems,cmdl.nodes(:,1),cmdl.nodes(:,2));
set(h,'Color',[0,0,1],'LineWidth',2);
hold off
axis image

print -dpng -r125 square_mesh02a.png
