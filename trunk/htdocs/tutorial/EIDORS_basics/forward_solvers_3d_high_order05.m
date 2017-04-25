%Difference between quadratic and linear approximation internal voltage
J0 = calc_jacobian( img0 );
J1 = calc_jacobian( img1 );
J2 = calc_jacobian( img2 );

for i=1:6;
   imgout = img0;
   switch i; case 1; imgout.elem_data = J0(4,:);
             case 2; imgout.elem_data = J1(4,:);
             case 3; imgout.elem_data = J2(4,:);
             case 4; continue
             case 5; imgout.elem_data = J0(4,:) - J1(4,:);
             case 6; imgout.elem_data = J0(4,:) - J2(4,:);
   end
   subplot(2,3,i);
   show_slices(imgout,[inf,inf,2.5]);
   eidors_colourbar(imgout);
end

print_convert forward_solvers_3d_high_order05a.png
