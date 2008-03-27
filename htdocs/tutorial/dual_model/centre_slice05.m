% 2D solver $Id: centre_slice05.m,v 1.1 2008-03-27 20:56:00 aadler Exp $

% Create and show inverse solver
imdl = mk_common_model('b3cr',[16,2]);
f_mdl= imdl.fwd_model;

% Create coarse model
imdl2d= mk_common_model('b2c2',16);
c_mdl= imdl2d.fwd_model;

% Show fine model
show_fem(f_mdl);
crop_model(gca, inline('x-z<-.5','x','y','z'))
view(-23,10)

scl= 1; % scale difference between c_mdl and f_mdl
c_els= c_mdl.elems;
c_ndsx= c_mdl.nodes(:,1)*scl;
c_ndsy= c_mdl.nodes(:,2)*scl;
c_ndsz= 0*c_ndsx; 

layersep= .3;
layerhig= .1;
hold on
% Lower resonstruction layer
hh= trimesh(c_els, c_ndsx, c_ndsy, c_ndsz+(-layersep)*scl);
set(hh, 'EdgeColor', [1,0,0]);
hh= trimesh(c_els, c_ndsx, c_ndsy, c_ndsz+(-layersep-layerhig)*scl);
set(hh, 'EdgeColor', [.3,.3,.3]);
hh= trimesh(c_els, c_ndsx, c_ndsy, c_ndsz+(-layersep+layerhig)*scl);
set(hh, 'EdgeColor', [.3,.3,.3]);

% Upper resonstruction layer
hh= trimesh(c_els, c_ndsx, c_ndsy, c_ndsz+(+layersep)*scl);
set(hh, 'EdgeColor', [0,0,1]);
hh= trimesh(c_els, c_ndsx, c_ndsy, c_ndsz+(+layersep-layerhig)*scl);
set(hh, 'EdgeColor', [.3,.3,.3]);
hh= trimesh(c_els, c_ndsx, c_ndsy, c_ndsz+(+layersep+layerhig)*scl);
set(hh, 'EdgeColor', [.3,.3,.3]);

[xs,ys,zs]=sphere(10);
for i=1:size(xyzr_pt,2)
   xp=xyzr_pt(1,i)/15;   yp=xyzr_pt(2,i)/15;
   zp=xyzr_pt(3,i)/15-1; rp=xyzr_pt(4,i)/15;
   hh=surf(rp*xs+xp, rp*ys+yp, rp*zs+zp);
   set(hh,'EdgeColor',[.4,0,.4]);
end

hold off;

print -r100 -dpng centre_slice05a.png;
