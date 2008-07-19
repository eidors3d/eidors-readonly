% Select ROI's

img = i_injury; img.elem_data = img.elem_data(:,630);
rimg = calc_colours( calc_slices( img ), img);

np= calc_colours('npoints');
xlocn= 0.25 * np;
ylocn= [2:5]/8 * np; 
for yl = 1:4;
   rimg(ylocn(yl) + (-2:2), xlocn + (-2:2) ) = 1;
end
image(rimg); axis square

for yl = 1:4;
   text(xlocn-1, ylocn(yl), num2str(yl));
end

axis off
print -r100 -dpng if_peep_trial03.png
