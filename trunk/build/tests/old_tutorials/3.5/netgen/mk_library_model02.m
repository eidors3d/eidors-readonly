fmdl= mk_library_model('adult_male_16el_lungs');
img = mk_image(fmdl, 0.25); % background conductivity
img.elem_data([fmdl.mat_idx{2};fmdl.mat_idx{3}]) = 0.05; % lungs
subplot(211); show_fem(img);
print_convert mk_library_model02a.png
view(2);
print_convert mk_library_model02b.png
