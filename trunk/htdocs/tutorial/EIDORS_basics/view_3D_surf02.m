imgv = rmfield(img,'elem_data');
imgv.node_data = vh.volt(:,1);

colours = calc_colours(imgv,[]);
patch('Faces',fmdl.boundary,'Vertices',fmdl.nodes, 'facecolor','interp', ...
      'facevertexcdata',colours,'CDataMapping','direct'); 

print_convert view_3D_surf02a.jpg '-density 75'
view(0,0)
print_convert view_3D_surf02b.jpg '-density 75'
