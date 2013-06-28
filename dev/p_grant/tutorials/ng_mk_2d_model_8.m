a = shape_library('get','pig_23kg','boundary');

% Fourier desriptor
pp = fourier_fit(a,length(a)-1); % avoid overfitting

% smoothened contour
ls = linspace(0,1,46); ls(end) = [];
a = fourier_fit(pp,ls);

% counter-clockwise contour
a = flipud(a);

% lungs, can also use ‘lungs_and_heart’ or ‘heart’
b = shape_library('get','pig_23kg','lungs');
pp = fourier_fit(b,length(b)-1); 
ls = linspace(0,1,46); ls(end) = [];
b = fourier_fit(pp,ls);
b = flipud(b);

% model with a hole
mdl1 = ng_mk_2d_model({a,b,0.05},[16,0.5],0.1);

%change netgen options to discourage additional nodes on contours
ng_write_opt('meshoptions.fineness',1,'options.meshsize',0.05);

% make second model of just lungs
mdl2 = ng_mk_2d_model(b);

% merge the two models
mdl = merge_meshes(mdl1,mdl2);

% make an image
img = mk_image(mdl,1);
img.elem_data(mdl.mat_idx{2}) = 0.3; 
show_fem(img,[0 1 0])

% remove the ng.opt file
delete('ng.opt');
print_convert eight.png