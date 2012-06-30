img.fwd_solve.get_all_meas = 1;
img.elem_data(fmdl.mat_idx{2}) = conduct; % Homogenous
vh = fwd_solve(img);

img.elem_data(fmdl.mat_idx{2}) = 0.95*conduct; %Non-conductive inclusion
vi = fwd_solve(img);

plot([vh.meas, 100*(vh.meas-vi.meas)])
axis tight
print_convert forward_solver_3d_03a.png '-density 77'
