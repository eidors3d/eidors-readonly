fmdl= mk_library_model('adult_male_16el_lungs');
img = mk_image(fmdl, 1); % background conductivity
img.elem_data(fmdl.mat_idx{2}) = 0.3001; % lungs
img.elem_data(fmdl.mat_idx{3}) = 0.3002; % lungs
ROI = calc_slices(img,[inf,inf,0.5]);
llung_ROI = ~isnan(ROI) & (ROI==0.3001);
rlung_ROI = ~isnan(ROI) & (ROI==0.3002);
thorax_ROI= ~isnan(ROI); % include lungs, too

subplot(131); imagesc(thorax_ROI); axis image
subplot(132); imagesc(rlung_ROI);  axis image
subplot(133); imagesc(llung_ROI);  axis image
print_convert('GREIT_IBEX_01a.png');

