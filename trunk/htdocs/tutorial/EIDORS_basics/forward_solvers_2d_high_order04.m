%Difference between quadratic and linear approximation internal voltage
img12n=img1n; 
img12n.node_data=v2all(1:size(fmdl.nodes,1),1)-v1all(1:size(fmdl.nodes,1),1);

%Plot the difference 
clf; subplot(221); show_fem(img12n,1);

%Difference between quadratic and linear approximation internal voltage
img13n=img1n; 
img13n.node_data=v3all(1:size(fmdl.nodes,1),1)-v1all(1:size(fmdl.nodes,1),1);

%Plot the difference 
subplot(222); show_fem(img13n,1);

print_convert forward_solvers_2d_high_order04a.png
