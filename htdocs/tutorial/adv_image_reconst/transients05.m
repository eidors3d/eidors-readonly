rimg = calc_slices(imgn);

slice = rimg(:,32,:);
image(calc_colours(reshape(slice,64,[]), imgn))

print_convert transients05a.png

