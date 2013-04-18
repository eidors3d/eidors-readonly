% get contours
thorax = shape_library('get','adult_male','boundary');
rlung  = shape_library('get','adult_male','right_lung');
llung  = shape_library('get','adult_male','left_lung');
% one could also run:
% shape_library('get','adult_male');
% to get all the info at once in a struct

% show the library image
shape_library('show','adult_male');
print_convert thoraxmdl01a.jpg '-density 100'
