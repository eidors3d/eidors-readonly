%Difference between quadratic and linear approximation internal voltage
J0 = calc_jacobian( img0 );
J1 = calc_jacobian( img1 );
%J2 = calc_jacobian( img2 ); -- curently doesn't work
%J3 = calc_jacobian( img3 ); 
img0j  = img0; img0j.elem_data = J0(4,:);
img01j = img0; img01j.elem_data = J1(4,:) - J0(4,:);

%Plot the difference 
clf;
subplot(221); show_fem(img0j,1);
subplot(223); show_fem(img0j,1); axis([-0.2, 0.6, 0.3, 1.1]);


%Plot the difference 
subplot(222); show_fem(img01j,1);
subplot(224); show_fem(img01j,1); axis([-0.2, 0.6, 0.3, 1.1]);

print_convert forward_solvers_2d_high_order05a.png
