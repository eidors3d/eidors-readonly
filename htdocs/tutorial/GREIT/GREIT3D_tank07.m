img = inv_solve(imdl3,vh,[vi{:}]);

img.show_slices.img_cols= 9;
levels = [inf,inf,1.5;
          inf,inf,2.0];
show_slices(img,levels);

print_convert GREIT3D_tank07a.jpg

img.elem_data = img.elem_data(:,4); % choose 4th image
show_fem(img);
print_convert GREIT3D_tank07b.jpg

show_3d_slices(img,[1.5,2],0,0); view(-16,18);
print_convert GREIT3D_tank07c.jpg
