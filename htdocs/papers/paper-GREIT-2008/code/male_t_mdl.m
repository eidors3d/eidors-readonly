function male_t_mdl( file );
% Build neonate thorax model from A. Tizzard's mesh
% $Id$
load HMT027.mat
fmdl.name= 'Thorax FEM';
fmdl.type= 'fwd_model';
fmdl.elems= tri;
fmdl.nodes= vtx;
fmdl.boundary= double( find_boundary(tri) );

%male
sep=[ ...
 -2.8560, -2.3340, -1.8887, -1.4587, -1.1516, -0.8906, -0.5681, -0.1996, ...
  0.1843,  0.5681,  0.8599,  1.1823,  1.5355,  1.9040,  2.3186,  2.7946]


eangl= atan2(vtx(enodes,1),vtx(enodes,2));
for i=1:16
   if i==1; 
      this_elec= find(eangl <  sep(1) | eangl >= sep(16)); 
   else
      this_elec= find(eangl >= sep(i-1) & eangl < sep(i)); 
   end
   fmdl.electrode(i).nodes = enodes(this_elec);
   fmdl.electrode(i).z_contact= 0.1;
end
fmdl.stimulation = [];
fmdl.gnd_node= 1;
show_fem(fmdl);
