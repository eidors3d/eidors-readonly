% $Id: square_mesh01.m,v 1.3 2008-04-17 19:45:28 aadler Exp $

z_contact= 0.01;
n_elec= 17;
nodes_per_elec= 5;
elec_width= 0.2;
elec_spacing= 1.0;

xllim=-12; xrlim= 12; ydepth=-15;
[x,y] = meshgrid( linspace(xllim,xrlim,49), linspace(ydepth,0,31) );
vtx= [x(:),y(:)];
% Refine points close to electrodes - don't worry if points overlap
[x,y] = meshgrid( -9:.25:9, -3:.25:0 );
vtx= [vtx; x(:),y(:)];

xgrid=  linspace(-elec_width/2, +elec_width/2, nodes_per_elec)';
x2grid= elec_width* [-5,-4,-3,3,4,5]'/4;
for i=1:n_elec
% Electrode centre
  x0= (i-1-(n_elec-1)/2)*elec_spacing;
  y0=0;
  elec_nodes{i}= [x0+ xgrid, y0+0*xgrid];
  vtx= [ vtx; ...
        [x0 + x2grid   ,  y0             + 0*x2grid];
        [x0 + xgrid*1.5,  y0-elec_width/2+ 0*xgrid];
        [x0 + x2grid*1.5, y0-elec_width/2+ 0*x2grid];
        [x0 + xgrid*2   , y0-elec_width  + 0*xgrid];
        [x0 + xgrid*2   , y0-elec_width*2+ 0*xgrid]];
end

fmdl= mk_fmdl_from_nodes( vtx, elec_nodes, z_contact, 'sq_m1');


subplot(121)
show_fem(fmdl); axis image
subplot(122)
show_fem(fmdl); axis image; axis([-2 2 -2.5 0.5]);

print -dpng -r150 square_mesh01a.png
