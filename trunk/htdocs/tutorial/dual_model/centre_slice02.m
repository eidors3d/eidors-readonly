% 2D solver $Id: centre_slice02.m,v 1.5 2008-03-28 02:28:30 aadler Exp $

% Create and show inverse solver
imdl = mk_common_model('b3cr',[16,2]);

load ng_mdl_16x2_coarse; f_mdl = ng_mdl_16x2_coarse;

imdl.fwd_model = f_mdl;

% Create coarse model
imdl2d= mk_common_model('b2c2',16);
c_mdl= imdl2d.fwd_model;

% Show fine model
show_fem(f_mdl);
crop_model(gca, inline('x-z<-15','x','y','z'))
view(-23,10)

scl= 15; % scale difference between c_mdl and f_mdl
c_els= c_mdl.elems;
c_ndsx= c_mdl.nodes(:,1)*scl;
c_ndsy= c_mdl.nodes(:,2)*scl;
c_ndsz= 0*c_ndsx; 

layersep= .3;
layerhig= .1;
hold on
% Lower resonstruction layer
hh= trimesh(c_els, c_ndsx, c_ndsy, c_ndsz+(1-layersep)*scl);
set(hh, 'EdgeColor', [1,0,0]);
hh= trimesh(c_els, c_ndsx, c_ndsy, c_ndsz+(1-layersep-layerhig)*scl);
set(hh, 'EdgeColor', [.3,.3,.3]);
hh= trimesh(c_els, c_ndsx, c_ndsy, c_ndsz+(1-layersep+layerhig)*scl);
set(hh, 'EdgeColor', [.3,.3,.3]);

% Upper resonstruction layer
hh= trimesh(c_els, c_ndsx, c_ndsy, c_ndsz+(1+layersep)*scl);
set(hh, 'EdgeColor', [0,0,1]);
hh= trimesh(c_els, c_ndsx, c_ndsy, c_ndsz+(1+layersep-layerhig)*scl);
set(hh, 'EdgeColor', [.3,.3,.3]);
hh= trimesh(c_els, c_ndsx, c_ndsy, c_ndsz+(1+layersep+layerhig)*scl);
set(hh, 'EdgeColor', [.3,.3,.3]);

for i=1:n_sims
   xp=xyzr_pt(1,i); yp=xyzr_pt(2,i);
   zp=xyzr_pt(3,i); rp=xyzr_pt(4,i);
   hh=surf(rp*xs+xp, rp*ys+yp, rp*zs+zp);
   set(hh,'EdgeColor',[.4,0,.4]);
end

hold off;

print -r100 -dpng centre_slice02a.png;
