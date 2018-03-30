   clear opt;
   opt.imgsz = [32 32];
   opt.square_pixels = true;
   opt.noise_figure = 0.5;
   img = mk_image(fmdl,1);
   imdl2b= mk_GREIT_model(img, 0.25, [], opt);

load mar28-horse_breathing.mat

for i=1:3; switch i;
   case 1; imdl = imdl3;  vv= horse3d; % 2 planes, 3D GREIT
   case 2; imdl = imdl2a; vv= horse2d; % 1 plane, 3D GREIT
   case 3; imdl = imdl2b; vv= horse2d; % 1 plane, 2D GREIT
   end
  
   subplot(3,2,2*i-1);
   img = inv_solve(imdl,vv(:,1),vv(:,2:end));
   img.calc_colours.ref_level = 0;
   show_slices(img,[inf,inf,2]);

   subplot(3,2,2*i-0);
   img.elem_data = img.elem_data(:,6);
   if i<3
   show_3d_slices(img,[1.6,2.0,2.4],[],[0.4]); view(-20,20);
   else
   show_slices(img); view(-20,20);
   end
end

print_convert mar29-horse-breathing.jpg

