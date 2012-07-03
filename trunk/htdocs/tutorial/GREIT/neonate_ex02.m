% solve with reference to the mean
imgall = inv_solve(imdl,mean(vv,2),vv);
[ins, exp] = find_frc(imgall,[],13,[],2); % find breaths


% use expirations as reference
vh = mean(vv(:,exp),2);  % reference is average
img = inv_solve(imdl,vh,vv);
img.elem_data = img.elem_data(:,ins); %only show inspirations
calc_colours('defaults');
calc_colours('ref_level',0);
calc_colours('backgnd',[1 1 1]);
calc_colours('greylev', 0.001);
img.show_slices.img_cols = 4;
img.show_slices.sep      = 2;
img.calc_colours.ref_level=0;
show_slices(img);

print_convert neonate_ex02a.png
