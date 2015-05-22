%% Simulation
% homogeneous
imgh = mk_image(fmdl,1);
vh = fwd_solve(imgh);
% spherical target
select_fun = inline('(x-.5).^2+(y-.5).^2+(z-1.25).^2<=0.1^2','x','y','z');
% inhomogeneous
imgi = imgh;
imgi.elem_data = imgh.elem_data + elem_select(fmdl,select_fun);
vi = fwd_solve(imgi);
show_fem(imgi)
print_convert vox_GREIT_sim_02a.png

clf
show_3d_slices(imgi,1.25,.5,.5);% cuts through the target center
print_convert vox_GREIT_sim_02b.png

clf
levels(:,3) = 2.25:-.5:.25; % cuts through the target center
levels(:,1:2) = Inf;
show_slices(imgi,levels)
print_convert vox_GREIT_sim_02c.png