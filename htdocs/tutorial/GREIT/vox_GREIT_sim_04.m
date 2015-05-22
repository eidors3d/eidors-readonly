show_fem(rimg)
print_convert vox_GREIT_sim_04a.png

clf
show_3d_slices(rimg,1.25,.5,.5);% cuts through the target center
print_convert vox_GREIT_sim_04b.png

clf
levels(:,3) = 2.5:-.5:.5; % cuts through the middle of the voxels
levels(:,1:2) = Inf;
show_slices(rimg,levels)
print_convert vox_GREIT_sim_04c.png