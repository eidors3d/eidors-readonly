%Get internal linear voltage distribution for first stimulation
v1all = v1.volt; 
img1n = rmfield(img1,'elem_data');
img1n.node_data = v1all(1:size(fmdl.nodes,1),1); %add first stim data

%Plot the distribution
subplot(121); show_slices(img1n,[inf,inf,2.5]);
eidors_colourbar(img1n);

%Get internal voltage distribution for difference eidors/high order
img11n=img1n; 
img11n.node_data=v1all(1:size(fmdl.nodes,1),1)-v0all(1:size(fmdl.nodes,1),1);

%Plot the difference of two linear approximations
subplot(122); show_slices(img11n,[inf,inf,2.5]);
eidors_colourbar(img11n);

print_convert forward_solvers_3d_high_order02a.png
