%Difference between quadratic and linear approximation internal voltage
img12n=img1n; 
img12n.node_data=v2all(1:size(fmdl.nodes,1),1)-v1all(1:size(fmdl.nodes,1),1);

subplot(121); show_slices(img1n,[inf,inf,2.5]);
eidors_colourbar(img1n);
%Plot the difference 
subplot(122); show_slices(img12n,[inf,inf,2.5]);
eidors_colourbar(img12n);

print_convert forward_solvers_3d_high_order04a.png
