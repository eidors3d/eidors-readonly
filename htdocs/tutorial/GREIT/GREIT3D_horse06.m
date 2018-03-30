load horse-breathing3D2D.mat

for i=1:3; switch i;
   case 1; imdl = imdl3;  vv= horse3d; % 2 planes, 3D GREIT
   case 2; imdl = imdl2a; vv= horse2d; % 1 plane, 3D GREIT
   case 3; imdl = imdl2b; vv= horse2d; % 1 plane, 2D GREIT
   end
  
   img = inv_solve(imdl,vv(:,1),vv(:,2:end));
   img.calc_colours.ref_level = 0;
   subplot(121);
   show_slices(img,[inf,inf,2]);

   subplot(122);
   img.elem_data = img.elem_data(:,6);
   if i<3; show_3d_slices(img,[1.6,2.0,2.4],[],[0.4]);
   else;   show_slices(img);
   end;    view(-20,20);

   print_convert(sprintf('GREIT3D_horse06%c.jpg',64+i));
end

