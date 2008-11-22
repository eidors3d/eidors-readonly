fmdl.name= 'Male Thorax FEM';
fmdl.type= 'fwd_model';
fmdl.elems= tri;
fmdl.nodes= vtx;
fmdl.boundary= find_boundary(tri);
for i= 1:size(fmdl.boundary,1)
   fmbi= fmdl.boundary(i,:);
   if any(fmbi(1) == enodes) &&  ...
      any(fmbi(2) == enodes) &&  ...
      any(fmbi(3) == enodes); 
      disp(i)
    end
end
fmdl.electrode = [];
%fmdl.electrode(1).nodes = enodes(1:8);
%fmdl.electrode(:).z_contact= 0.1;
fmdl.stimulation = [];
fmdl.gnd_node= 1;
show_fem(fmdl);
for i=101:117;xyz=vtx(enodes(i),:);line(xyz(1),xyz(2),xyz(3),'marker','o','markerfacecolor','cyan','markersize',10); pause; end

