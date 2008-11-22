% Build neonate thorax model from A. Tizzard's mesh
% $Id$
clear fmdl;
fmdl.name= 'Thorax FEM';
fmdl.type= 'fwd_model';
fmdl.elems= tri;
fmdl.nodes= vtx;
fmdl.boundary= double( find_boundary(tri) );

% female
sep=[ ...
 -2.8560, -2.4415, -2.1036, -1.7505, -1.4587, -1.1209, -0.7370, -0.3378, ...
  0.2150,  0.6756,  1.0749,  1.4127,  1.7658,  2.1497,  2.5489,  2.9942];


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
%show_fem(fmdl);
save female_t_mdl.mat fmdl
