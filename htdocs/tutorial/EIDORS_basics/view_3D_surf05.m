imgv = rmfield(img,'elem_data');
imgv.node_data = vh.volt(:,1);

imgv.calc_colours.clim = 0.3; % Colour limits
colours = calc_colours(imgv,[]);
show_fem(fmdl);
patch('Faces',fmdl.boundary,'Vertices',fmdl.nodes, 'facecolor','interp', ...
      'facevertexcdata',colours,'CDataMapping','direct'); 

print_convert view_3D_surf05a.png '-density 75'
view(0,0)
print_convert view_3D_surf05b.png '-density 75'
