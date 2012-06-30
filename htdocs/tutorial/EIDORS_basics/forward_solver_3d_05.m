img.fwd_model.electrode([2,13]) = img.fwd_model.electrode([13,2]); % flip electrodes

img.elem_data(fmdl.mat_idx{2}) = conduct; % Homogenous
vh = fwd_solve(img);

img.elem_data(fmdl.mat_idx{2}) = 0.95*conduct; %Non-conductive inclusion
vi = fwd_solve(img);

img_v.node_data = vh.volt(:,1);

show_slices(img_v,[0.25;0.35;0.45]*[inf,inf,1])
print_convert forward_solver_3d_05a.png '-density 75'

img_v.node_data = vh.volt(:,1) - vi.volt(:,1);
show_slices(img_v,[0.25;0.35;0.45]*[inf,inf,1])
print_convert forward_solver_3d_05b.png '-density 75'
