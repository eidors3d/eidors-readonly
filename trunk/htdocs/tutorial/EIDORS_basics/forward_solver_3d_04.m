img_v = rmfield(img, 'elem_data');
img_v.node_data = vh.volt(:,1);

show_slices(img_v,[0.2;0.3;0.4]*[inf,inf,1])
print_convert forward_solver_3d_04a.png '-density 75'

img_v.node_data = vh.volt(:,1) - vi.volt(:,1);
show_slices(img_v,[0.2;0.3;0.4]*[inf,inf,1])
print_convert forward_solver_3d_04b.png '-density 75'
