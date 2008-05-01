% $Id: cube_mesh01.m,v 1.2 2008-05-01 16:27:44 aadler Exp $

z_contact= 0.01;
nx_elec= 4;
ny_elec= 4;
nodes_per_elec= 3;
elec_width= 0.2;
elec_spacing= 1.0;

xllim=-8; xrlim= 8;
yllim=-8; yrlim= 8; zdepth=-10;
[x,y,z] = meshgrid( linspace(xllim,xrlim,32+1), ...
                    linspace(yllim,yrlim,32+1), ...
                    linspace(zdepth,0,20+1) );
vtx= [x(:),y(:),z(:)];
% Refine points close to electrodes - don't worry if points overlap
[x,y,z] = meshgrid( -4:.25:4, -4:.25:4, -2:.25:0 );
vtx= [vtx; x(:),y(:),z(:)];

nx_elec= 1;
ny_elec= 1;
nodes_per_elec= 2;
xllim=-2; xrlim= 2;
yllim=-2; yrlim= 2; zdepth=-1;
[x,y,z] = meshgrid( linspace(xllim,xrlim,8+1), ...
                    linspace(yllim,yrlim,8+1), ...
                    linspace(zdepth,0,2+1) );
vtx= [x(:),y(:),z(:)];

% Define each electrode
xy_elec_grid= linspace( -elec_width/2,elec_width/2, nodes_per_elec);
[xe,ye,ze] = meshgrid( xy_elec_grid, xy_elec_grid, 0);
xy_elec_grid= linspace( -elec_width,elec_width, 2*nodes_per_elec-1);
[x,y,z] = meshgrid( xy_elec_grid, xy_elec_grid, -1:.25:0);

k=0;
for i= -(nx_elec-1)/2:(nx_elec-1)/2
  for j= -(ny_elec-1)/2:(ny_elec-1)/2
% Electrode centre
     x0= i*elec_spacing; y0=j*elec_spacing; z0=0;
     k=k+1;elec_nodes{k}= [x0+ xe(:), y0+ye(:), z0+ze(:)];
%    vtx= [ vtx; [x0+x(:), y0+y(:), z0+z(:)]];
   end
end

%vtx= vtx+ .00001*randn(size(vtx));
fmdl= mk_fmdl_from_nodes( vtx, elec_nodes, z_contact, 'sq_m1');


%subplot(121)
show_fem(fmdl); axis image
%subplot(122)
%show_fem(fmdl); axis image; axis([-2 2 -2.5 0.5]);

print -dpng -r150 cube_mesh01a.png
