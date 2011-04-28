fmdl= mk_library_model('adult_male_16el_lungs');
fmdl.electrode = fmdl.electrode([9:16,1:8]);
img = mk_image(fmdl, 1); % background conductivity
% Note, I don't understand why it takes a conductivity increase
%  here: This needs to be debugged further
img.elem_data([fmdl.mat_idx{2};fmdl.mat_idx{3}]) = 1.3; % lungs
subplot(211); show_fem(img);
print_convert adult_ex01a.png
show_fem(img,[0,1]); view(2); %electrode #'s
print_convert adult_ex01b.png
