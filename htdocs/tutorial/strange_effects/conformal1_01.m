xllim=-12; xrlim= 12; ydepth=-15;
[x,y] = meshgrid( linspace(xllim,xrlim,49), linspace(ydepth,0,31) );
vtx= [x(:),y(:)];

elec_nodes{1}= [x(1,:);y(1,:)]';
elec_nodes{2}= [x(end,:);y(end,:)]';

z_contact= 0.01;
fmdl= mk_fmdl_from_nodes( vtx, elec_nodes, z_contact, 'sq_m1');

subplot(121)
show_fem(fmdl); axis image

