% Geophysics model $Id: dm_geophys01.m,v 1.2 2008-04-03 19:33:27 aadler Exp $

n_nodes= 1000;
z_contact= 0.01;
n_elec= 9;
nodes_per_elec= 5;
elec_width= 0.2;
elec_spacing= 1.0;

xllim=-10; xrlim= 10; ydepth=-10;
xgrid=  linspace(-elec_width/2, +elec_width/2, nodes_per_elec)';
x2grid= elec_width* [-5,-4,-3,3,4,5]'/4;
refine_nodes= [xllim,0; xrlim,0; 0, ydepth];
for i=1:n_elec
% Electrode centre
  x0= (i-1-(n_elec-1)/2)*elec_spacing;
  y0=0;
  elec_nodes{i}= [x0+ xgrid, y0+0*xgrid];
  refine_nodes= [ refine_nodes; ...
                 [x0 + x2grid   ,  y0             + 0*x2grid];
                 [x0 + xgrid*1.5,  y0-elec_width/2+ 0*xgrid];
                 [x0 + x2grid*1.5, y0-elec_width/2+ 0*x2grid];
                 [x0 + xgrid*2   , y0-elec_width  + 0*xgrid];
                 [x0 + xgrid*2   , y0-elec_width*2+ 0*xgrid]];
end
  
% Distance to edge
%fd=inline('-min([10+p(:,1),10-p(:,1),-p(:,2),10+p(:,2)],[],2)','p');
fd=inline('-min(-p(:,2),10-sqrt(sum(p.^2,2)))','p');
bbox = [xllim,ydepth;xrlim,0];
gmdl= dm_mk_fwd_model( fd, [], n_nodes, bbox, ...
                          elec_nodes, refine_nodes, z_contact);


% Refine elements near upper boundary
n_nodes= 3000;
fh=inline('min(0.5-p(:,2)/4,2)','p');
fmdl= dm_mk_fwd_model( fd, fh, n_nodes, bbox, ...
                          elec_nodes, refine_nodes, z_contact);

% Show results - big
subplot(121)
show_fem(gmdl); axis on
subplot(122)
show_fem(fmdl); axis on

print -dpng -r150 dm_geophys01a.png

% Show results - small
subplot(121)
show_fem(gmdl); axis([-3,3,-3,0.3]);
subplot(122)
show_fem(fmdl); axis([-3,3,-3,0.3]);
print -dpng -r150 dm_geophys01b.png
