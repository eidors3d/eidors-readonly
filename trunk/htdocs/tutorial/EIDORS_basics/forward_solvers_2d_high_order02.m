%Get internal linear voltage distribution for first stimulation
v1all = v1.volt; 
img1n = rmfield(img1,'elem_data');
img1n.node_data = v1all(1:size(fmdl.nodes,1),1); %add first stim data

%Plot the distribution
img1n.calc_colours.cb_shrink_move = [0.5,0.5,0];
clf; subplot(221); show_fem(img1n,1);

%Get internal voltage distribution for difference eidors/high order
img11n=img1n; 
img11n.node_data=v1all(1:size(fmdl.nodes,1),1)-v0all(1:size(fmdl.nodes,1),1);

%Plot the difference of two linear approximations
img11n.calc_colours.cb_shrink_move = [0.5,0.5,0];
subplot(222); show_fem(img11n,1);

print_convert forward_solvers_2d_high_order02a.png
